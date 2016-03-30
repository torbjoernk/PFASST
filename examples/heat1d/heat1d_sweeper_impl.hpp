#include "heat1d_sweeper.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include <vector>
using std::shared_ptr;
using std::vector;

#include <pfasst/globals.hpp>
#include <pfasst/util.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>


namespace pfasst
{
  namespace examples
  {
    namespace heat1d
    {
      template<class SweeperTrait, typename Enabled>
      void
      Heat1D<SweeperTrait, Enabled>::init_opts()
      {
        config::options::add_option<size_t>("Heat 1D", "num_dofs",
                                            "number spatial degrees of freedom on fine level");
        config::options::add_option<size_t>("Heat 1D", "coarse_factor",
                                            "coarsening factor");
        config::options::add_option<spatial_t>("Heat 1D", "nu",
                                               "thermal diffusivity");
      }

      template<class SweeperTrait, typename Enabled>
      Heat1D<SweeperTrait, Enabled>::Heat1D(const size_t ndofs)
        :   IMEX<SweeperTrait, Enabled>()
          , _lap(ndofs)
      {
        this->encap_factory()->set_size(ndofs);

        for (size_t i = 0; i < ndofs; ++i) {
          spatial_t kx = TWO_PI
                            * ((i <= ndofs / 2) ? spatial_t(i)
                                                : spatial_t(i) - spatial_t(ndofs));
          this->_lap[i] = pfasst::almost_zero(pow(kx, 2)) ? 0.0 : -pow(kx, 2);
        }
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat1D<SweeperTrait, Enabled>::set_options()
      {
        IMEX<SweeperTrait, Enabled>::set_options();

        this->_nu = config::get_value<spatial_t>("nu", this->_nu);
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat1D<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();

        // taken from pySDC
        //   xvalues = np.array([(i) * self.dx for i in range(self.nvars)])
        //   me.values = np.sin(2 * np.pi * xvalues) * np.exp(-t * (2 * np.pi)**2 * self.nu)
        const spatial_t dx = 1.0 / spatial_t(this->get_num_dofs());
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          result->data()[i] = std::sin(TWO_PI * i * dx) * exp(-t * 4 * PI_SQR * this->_nu);
        }

        return result;
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat1D<SweeperTrait, Enabled>::post_step()
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
      Heat1D<SweeperTrait, Enabled>::converged(const bool pre_check)
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
                   "Observables after "
                   << ((this->get_status()->get_iteration() == 0)
                          ? std::string("prediction")
                          : std::string("iteration ") + std::to_string(this->get_status()->get_iteration())));
          for (size_t m = 0; m < num_nodes; ++m) {
            ML_CVLOG(1, this->get_logger_id(),
                     "  t["<<m<<"]=" << LOG_FIXED << (t + dt * nodes[m])
                     << "      |abs residual| = " << LOG_FLOAT << this->_abs_res_norms[m]
                     << "      |rel residual| = " << LOG_FLOAT << this->_rel_res_norms[m]
                     << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[m])
                     << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[m]));
          }
          ML_CLOG(INFO, this->get_logger_id(),
                  "  t["<<num_nodes<<"]=" << LOG_FIXED << (t + dt * nodes[num_nodes])
                  << "      |abs residual| = " << LOG_FLOAT << this->_abs_res_norms[num_nodes]
                  << "      |rel residual| = " << LOG_FLOAT << this->_rel_res_norms[num_nodes]
                  << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[num_nodes])
                  << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[num_nodes]));
        }

        return converged;
      }

      template<class SweeperTrait, typename Enabled>
      bool
      Heat1D<SweeperTrait, Enabled>::converged()
      {
        return this->converged(false);
      }

      template<class SweeperTrait, typename Enabled>
      size_t
      Heat1D<SweeperTrait, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory().size();
      }


      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat1D<SweeperTrait, Enabled>::compute_error(const typename SweeperTrait::time_t& t)
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
          ML_CVLOG(3, this->get_logger_id(), LOG_FIXED << "error t=" << (t + ds) << ": "
                                          << LOG_FLOAT << to_string(error[m]));
        }

        return error;
      }

      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat1D<SweeperTrait, Enabled>::compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
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

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat1D<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        UNUSED(u);
        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "evaluating EXPLICIT part at t=" << t);

        auto result = this->get_encap_factory().create();
        result->zero();

        this->_num_expl_f_evals++;

        return result;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat1D<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "evaluating IMPLICIT part at t=" << t);

        spatial_t c = this->_nu / spatial_t(this->get_num_dofs());

        auto* z = this->_fft.forward(u);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] *= c * this->_lap[i];
        }

        auto result = this->get_encap_factory().create();
        this->_fft.backward(result);

        this->_num_impl_f_evals++;

        return result;
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat1D<SweeperTrait, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs)
      {
        ML_CVLOG(4, this->get_logger_id(),
                 LOG_FIXED << "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);

        spatial_t c = this->_nu * dt;

        auto* z = this->_fft.forward(rhs);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] /= (1.0 - c * this->_lap[i]) * spatial_t(this->get_num_dofs());
        }
        this->_fft.backward(u);

        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          f->data()[i] = (u->get_data()[i] - rhs->get_data()[i]) / dt;
        }

        this->_num_impl_solves++;
      }
    }  // ::pfasst::examples::heat1d
  }  // ::pfasst::examples
}  // ::pfasst
