#include "advec_diff_sweeper.hpp"

#include <cmath>
#include <complex>
using namespace std;

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
      template<class SweeperTrait, typename Enabled>
      void
      AdvecDiff<SweeperTrait, Enabled>::init_opts()
      {
        config::options::add_option<size_t>("Advection-Diffusion", "num_dofs", "number spatial degrees of freedom on fine level");
        config::options::add_option<size_t>("Advection-Diffusion", "coarse_factor", "coarsening factor");
        config::options::add_option<spatial_type>("Advection-Diffusion", "nu", "diffusivity");
        config::options::add_option<spatial_type>("Advection-Diffusion", "vel", "velocity of advection");
      }

      template<class SweeperTrait, typename Enabled>
      AdvecDiff<SweeperTrait, Enabled>::AdvecDiff(const size_t& ndofs, const typename SweeperTrait::spatial_type& nu, const typename SweeperTrait::spatial_type& v)
        :   IMEX<SweeperTrait, Enabled>()
          , _t0(1.0)
          , _nu(nu)
          , _v(v)
          , _ddx(ndofs)
          , _lap(ndofs)
      {
        this->encap_factory()->set_size(ndofs);

        for (size_t i = 0; i < ndofs; ++i) {
          spatial_type kx = two_pi<spatial_type>()
                            * ((i <= ndofs / 2) ? spatial_type(i)
                                                : spatial_type(i) - spatial_type(ndofs));
          this->_ddx[i] = complex<spatial_type>(0.0, 1.0) * kx;
          this->_lap[i] = pfasst::almost_zero(kx * kx) ? 0.0 : -kx * kx;
        }
      }

      template<class SweeperTrait, typename Enabled>
      void
      AdvecDiff<SweeperTrait, Enabled>::set_options()
      {
        IMEX<SweeperTrait, Enabled>::set_options();

        this->_nu = config::get_value<typename traits::spatial_type>("nu", DEFAULT_DIFFUSIVITY);
        this->_v = config::get_value<typename traits::spatial_type>("vel", DEFAULT_VELOCITY);
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      AdvecDiff<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_type& t)
      {
        auto result = this->get_encap_factory()->create();

        const spatial_type dx = 1.0 / sqrt(4.0 * pi<spatial_type>() * this->_nu * (t + this->_t0));
        
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          result->data()[i] = 0.0;
        }

        for (int ii = -2; ii < 3; ++ii) {
          for (size_t i = 0; i < this->get_num_dofs(); ++i) {
            spatial_type x = spatial_type(i) / this->get_num_dofs() - 0.5 + ii - t * this->_v;
            result->data()[i] += dx * exp(-x * x / (4 * this->_nu * (t + this->_t0)));
          }
        }

        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "EXACT t=" << t << ": " << LOG_FLOAT << to_string(result));

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

      template<class SweeperTrait, typename Enabled>
      bool
      AdvecDiff<SweeperTrait, Enabled>::converged()
      {
        const bool converged = IMEX<SweeperTrait, Enabled>::converged();

        assert(this->get_status() != nullptr);
        const time_type t = this->get_status()->get_time();
        const time_type dt = this->get_status()->get_dt();

        auto error = this->compute_error(t);
        auto rel_error = this->compute_relative_error(error, t);

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), time_type(0.0));

        ML_CVLOG(1, this->get_logger_id(), "Observables after " << ((this->get_status()->get_iteration() == 0) ? string("prediction") : string("iteration ") + to_string(this->get_status()->get_iteration())));
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

        return converged;
      }

      template<class SweeperTrait, typename Enabled>
      size_t
      AdvecDiff<SweeperTrait, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory()->size();
      }


      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_type>>
      AdvecDiff<SweeperTrait, Enabled>::compute_error(const typename SweeperTrait::time_type& t)
      {
        ML_CVLOG(4, this->get_logger_id(), "computing error");

        assert(this->get_status() != nullptr);
        const time_type dt = this->get_status()->get_dt();

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), time_type(0.0));

        vector<shared_ptr<encap_type>> error;
        error.resize(num_nodes + 1);
        generate(error.begin(), error.end(),
                 bind(&encap_type::factory_type::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          const time_type ds = dt * (nodes[m] - nodes[0]);
          error[m] = pfasst::encap::axpy(-1.0, this->exact(t + ds), this->get_states()[m]);
          ML_CVLOG(3, this->get_logger_id(), LOG_FIXED << "error t=" << (t + ds) << ": "
                                          << LOG_FLOAT << to_string(error[m]));
        }

        return error;
      }

      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_type>>
      AdvecDiff<SweeperTrait, Enabled>::compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_type>>& error,
                                                               const typename SweeperTrait::time_type& t)
      {
        UNUSED(t);

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), time_type(0.0));

        vector<shared_ptr<encap_type>> rel_error;
        rel_error.resize(error.size());
        generate(rel_error.begin(), rel_error.end(),
                 bind(&encap_type::factory_type::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          rel_error[m]->scaled_add(1.0 / this->get_states()[m]->norm0(), error[m]);
        }

        return rel_error;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      AdvecDiff<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_type& t,
                                                       const shared_ptr<typename SweeperTrait::encap_type> u)
      {
        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "evaluating EXPLICIT part at t=" << t);
        ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tu:   " << to_string(u));

        spatial_type c = - this->_v / spatial_type(this->get_num_dofs());

        auto* z = this->_fft.forward(u);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] *= c * this->_ddx[i];
        }

        auto result = this->get_encap_factory()->create();
        this->_fft.backward(result);

        this->_num_expl_f_evals++;

        ML_CVLOG(5, this->get_logger_id(), "\t  -> " << to_string(result));
        return result;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      AdvecDiff<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                       const shared_ptr<typename SweeperTrait::encap_type> u)
      {
        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "evaluating IMPLICIT part at t=" << t);
        ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tu:   " << to_string(u));

        spatial_type c = this->_nu / spatial_type(this->get_num_dofs());

        auto* z = this->_fft.forward(u);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] *= c * this->_lap[i];
        }

        auto result = this->get_encap_factory()->create();
        this->_fft.backward(result);

        this->_num_impl_f_evals++;

        ML_CVLOG(5, this->get_logger_id(), "\t  -> " << to_string(result));
        return result;
      }

      template<class SweeperTrait, typename Enabled>
      void
      AdvecDiff<SweeperTrait, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_type> f,
                                                    shared_ptr<typename SweeperTrait::encap_type> u,
                                                    const typename SweeperTrait::time_type& t,
                                                    const typename SweeperTrait::time_type& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_type> rhs)
      {
        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);
        ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tf:   " << to_string(f));
        ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tu:   " << to_string(u));
        ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\trhs: " << to_string(rhs));

        spatial_type c = this->_nu * dt;

        auto* z = this->_fft.forward(rhs);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] /= (1.0 - c * this->_lap[i]) * spatial_type(this->get_num_dofs());
        }
        this->_fft.backward(u);

        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          f->data()[i] = (u->get_data()[i] - rhs->get_data()[i]) / dt;
        }

        this->_num_impl_solves++;

        ML_CVLOG(5, this->get_logger_id(), "\t->");
        ML_CVLOG(5, this->get_logger_id(), "\t  f: " << to_string(f));
        ML_CVLOG(5, this->get_logger_id(), "\t  u: " << to_string(u));
      }
    }  // ::pfasst::examples::advec_diff
  }  // ::pfasst::examples
}  // ::pfasst
