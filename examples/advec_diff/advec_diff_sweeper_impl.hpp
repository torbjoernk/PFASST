#include "advec_diff_sweeper.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <memory>
#include <string>
#include <vector>
using std::shared_ptr;
using std::vector;

#include <leathers/push>
#include <leathers/all>
#include <boost/math/constants/constants.hpp>
#include <leathers/pop>
using boost::math::constants::pi;
using boost::math::constants::two_pi;
using boost::math::constants::pi_sqr;

#include <pfasst/globals.hpp>
#include <pfasst/util.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>


namespace pfasst
{
  namespace examples
  {
    namespace advec_diff
    {
      /**
       * | option          | description                                              |
       * |-----------------|----------------------------------------------------------|
       * | `num_dofs`      | number spatial degrees of freedom on finest level        |
       * | `coarse_factor` | coarsening factor from one level to the next coarser one |
       * | `nu`            | diffusivity coefficient                                  |
       * | `vel`           | velocity of advection                                    |
       */
      template<class SweeperTrait, typename Enabled>
      void
      AdvecDiff<SweeperTrait, Enabled>::init_opts()
      {
        config::options::add_option<size_t>("Advection-Diffusion", "num_dofs", "number spatial degrees of freedom on fine level");
        config::options::add_option<size_t>("Advection-Diffusion", "coarse_factor", "coarsening factor");
        config::options::add_option<typename traits::spatial_t>("Advection-Diffusion", "nu", "diffusivity");
        config::options::add_option<typename traits::spatial_t>("Advection-Diffusion", "vel", "velocity of advection");
      }

      template<class SweeperTrait, typename Enabled>
      AdvecDiff<SweeperTrait, Enabled>::AdvecDiff(const size_t& ndofs, const typename SweeperTrait::spatial_t& nu,
                                                  const typename SweeperTrait::spatial_t& v)
        :   IMEX<SweeperTrait, Enabled>()
          , _t0(1.0)
          , _nu(nu)
          , _v(v)
          , _ddx(ndofs)
          , _lap(ndofs)
      {
        this->encap_factory()->set_size(ndofs);

        for (size_t i = 0; i < ndofs; ++i) {
          typename traits::spatial_t kx = two_pi<typename traits::spatial_t>()
                            * ((i <= ndofs / 2) ? typename traits::spatial_t(i)
                                                : typename traits::spatial_t(i) - typename traits::spatial_t(ndofs));
          this->_ddx[i] = std::complex<typename traits::spatial_t>(0.0, 1.0) * kx;
          this->_lap[i] = pfasst::almost_zero(kx * kx) ? 0.0 : -kx * kx;
        }
      }

      template<class SweeperTrait, typename Enabled>
      void
      AdvecDiff<SweeperTrait, Enabled>::set_options()
      {
        IMEX<SweeperTrait, Enabled>::set_options();

        this->_nu = config::get_value<typename traits::spatial_t>("nu", DEFAULT_DIFFUSIVITY);
        this->_v = config::get_value<typename traits::spatial_t>("vel", DEFAULT_VELOCITY);
      }

      /**
       * The exact solution is computed on the spatial interval @f$ [-2.5, 2.5] @f$ as an
       * overlay of five Gaussian bells with the formula for each spatial component:
       * @f[
       *   u_i = \sum^{2}_{a=-2} \frac{1}{\sqrt{4 \pi \nu (t + t_0)}} e^{- \left(\frac{i}{N}
       *            - 0.5 + a - t * v \right)^2 / \left( 4 \nu (t + t_0) \right)}
       * @f]
       */
      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      AdvecDiff<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();

        const typename traits::spatial_t dx = 1.0 / sqrt(4.0 * pi<typename traits::spatial_t>() * this->_nu * (t + this->_t0));

        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          result->data()[i] = 0.0;
        }

        for (int ii = -2; ii < 3; ++ii) {
          for (size_t i = 0; i < this->get_num_dofs(); ++i) {
            typename traits::spatial_t x = typename traits::spatial_t(i) / this->get_num_dofs() - 0.5 + ii - t * this->_v;
            result->data()[i] += dx * exp(-x * x / (4 * this->_nu * (t + this->_t0)));
          }
        }

//         ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "EXACT t=" << t << ": " << LOG_FLOAT << to_string(result));

        return result;
      }

      template<class SweeperTrait, typename Enabled>
      void
      AdvecDiff<SweeperTrait, Enabled>::post_step()
      {
        IMEX<SweeperTrait, Enabled>::post_step();

        ML_CLOG(INFO, this->get_logger_id(), "number function evaluations:");
        ML_CLOG(INFO, this->get_logger_id(), "  expl:        " << this->_num_expl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl:        " << this->_num_impl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl solves: " << this->_num_impl_solves);

        this->_num_expl_f_evals = 0;
        this->_num_impl_f_evals = 0;
        this->_num_impl_solves = 0;
      }

      /**
       * If @p pre_check is `false`, it computes the absolute and relative errors at all time nodes
       * logs the norms of the residuals and errors at all time nodes.
       */
      template<class SweeperTrait, typename Enabled>
      bool
      AdvecDiff<SweeperTrait, Enabled>::converged(const bool pre_check)
      {
        const bool converged = IMEX<SweeperTrait, Enabled>::converged(pre_check);

        if (!pre_check) {
          assert(this->get_status() != nullptr);
          const typename traits::time_t t = this->get_status()->get_time();
          const typename traits::time_t dt = this->get_status()->get_dt();

          auto error = this->compute_error(t);
          auto rel_error = this->compute_relative_error(error, t);

          assert(this->get_quadrature() != nullptr);
          auto nodes = this->get_quadrature()->get_nodes();
          const auto num_nodes = this->get_quadrature()->get_num_nodes();
          nodes.insert(nodes.begin(), typename traits::time_t(0.0));

          ML_CVLOG(1, this->get_logger_id(),
                   "Observables after " << ((this->get_status()->get_iteration() == 0)
                                            ? std::string("prediction")
                                            : std::string("iteration ") + std::to_string(this->get_status()->get_iteration())));
          for (size_t m = 0; m < num_nodes; ++m) {
            ML_CVLOG(1, this->get_logger_id(), "  t["<<m<<"]=" << LOG_FIXED << (t + dt * nodes[m])
                                            << "      |abs residual| = " << LOG_FLOAT << this->_abs_res_norms[m]
                                            << "      |rel residual| = " << LOG_FLOAT << this->_rel_res_norms[m]
                                            << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[m])
                                            << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[m]));
          }
          ML_CLOG(INFO, this->get_logger_id(), "  t["<<num_nodes<<"]=" << LOG_FIXED << (t + dt * nodes[num_nodes])
                                            << "      |abs residual| = " << LOG_FLOAT << this->_abs_res_norms[num_nodes]
                                            << "      |rel residual| = " << LOG_FLOAT << this->_rel_res_norms[num_nodes]
                                            << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[num_nodes])
                                            << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[num_nodes]));
        }

        return converged;
      }

      template<class SweeperTrait, typename Enabled>
      bool
      AdvecDiff<SweeperTrait, Enabled>::converged()
      {
        return this->converged(false);
      }

      template<class SweeperTrait, typename Enabled>
      size_t
      AdvecDiff<SweeperTrait, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory().size();
      }


      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      AdvecDiff<SweeperTrait, Enabled>::compute_error(const typename SweeperTrait::time_t& t)
      {
        ML_CVLOG(4, this->get_logger_id(), "computing error");

        assert(this->get_status() != nullptr);
        const typename traits::time_t dt = this->get_status()->get_dt();

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), typename traits::time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> error;
        error.resize(num_nodes + 1);
        std::generate(error.begin(), error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          const typename traits::time_t ds = dt * (nodes[m] - nodes[0]);
          error[m] = pfasst::encap::axpy(-1.0, this->exact(t + ds), this->get_states()[m]);
//           ML_CVLOG(3, this->get_logger_id(), LOG_FIXED << "error t=" << (t + ds) << ": "
//                                           << LOG_FLOAT << to_string(error[m]));
        }

        return error;
      }

      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      AdvecDiff<SweeperTrait, Enabled>::compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                                               const typename SweeperTrait::time_t& t)
      {
        UNUSED(t);

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), typename traits::time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> rel_error;
        rel_error.resize(error.size());
        std::generate(rel_error.begin(), rel_error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          rel_error[m]->scaled_add(1.0 / this->get_states()[m]->norm0(), error[m]);
        }

        return rel_error;
      }

      /**
       * @details Evaluates the explicit part of right hand side of the problem equation - the
       *   advection - in Fourier space.
       *   I.e., the spatial data @p u is transformed into Fourier space and propagated in time by
       *    applying the Fourier discretization of the gradient w.r.t. to the velocity
       *    `AdvecDiff::_v`.
       */
      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      AdvecDiff<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "evaluating EXPLICIT part at t=" << t);
//         ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tu:   " << to_string(u));

        typename traits::spatial_t c = - this->_v / typename traits::spatial_t(this->get_num_dofs());

        auto* z = this->_fft.forward(u);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] *= c * this->_ddx[i];
        }

        auto result = this->get_encap_factory().create();
        this->_fft.backward(result);

        this->_num_expl_f_evals++;

//         ML_CVLOG(5, this->get_logger_id(), "\t  -> " << to_string(result));
        return result;
      }

      /**
       * @details Evaluates the implicit part of right hand side of the problem equation - the
       *   diffusion - in Fourier space.
       *   I.e., the spatial data @p u is transformed into Fourier space and propagated in time by
       *   applying the Fourier discretization of the Laplacian w.r.t. to the diffusivity
       *   coefficient `AdvecDiff::_nu`.
       */
      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      AdvecDiff<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "evaluating IMPLICIT part at t=" << t);
//         ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tu:   " << to_string(u));

        typename traits::spatial_t c = this->_nu / typename traits::spatial_t(this->get_num_dofs());

        auto* z = this->_fft.forward(u);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] *= c * this->_lap[i];
        }

        auto result = this->get_encap_factory().create();
        this->_fft.backward(result);

        this->_num_impl_f_evals++;

//         ML_CVLOG(5, this->get_logger_id(), "\t  -> " << to_string(result));
        return result;
      }

      /**
       * @details The implicit part here is the discretization of the Laplacian of the diffusion
       *    in the problem equation.
       *    This is solved in Fourier space.
       */
      template<class SweeperTrait, typename Enabled>
      void
      AdvecDiff<SweeperTrait, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs)
      {
        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);
//         ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tf:   " << to_string(f));
//         ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tu:   " << to_string(u));
//         ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\trhs: " << to_string(rhs));

        typename traits::spatial_t c = this->_nu * dt;

        auto* z = this->_fft.forward(rhs);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] /= (1.0 - c * this->_lap[i]) * typename traits::spatial_t(this->get_num_dofs());
        }
        this->_fft.backward(u);

        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          f->data()[i] = (u->get_data()[i] - rhs->get_data()[i]) / dt;
        }

        this->_num_impl_solves++;

//         ML_CVLOG(5, this->get_logger_id(), "\t->");
//         ML_CVLOG(5, this->get_logger_id(), "\t  f: " << to_string(f));
//         ML_CVLOG(5, this->get_logger_id(), "\t  u: " << to_string(u));
      }
    }  // ::pfasst::examples::advec_diff
  }  // ::pfasst::examples
}  // ::pfasst
