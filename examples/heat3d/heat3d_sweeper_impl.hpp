#include "heat3d_sweeper.hpp"

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
    namespace heat3d
    {
      template<class SweeperTrait, typename Enabled>
      void
      Heat3D<SweeperTrait, Enabled>::init_opts()
      {
        config::options::add_option<size_t>("Heat 3D", "num_dofs", "number spatial degrees of freedom per dimension on fine level");
        config::options::add_option<size_t>("Heat 3D", "coarse_factor", "coarsening factor");
        config::options::add_option<spatial_type>("Heat 3D", "nu", "thermal diffusivity");
      }

      template<class SweeperTrait, typename Enabled>
      Heat3D<SweeperTrait, Enabled>::Heat3D(const size_t& ndofs, const typename SweeperTrait::spatial_type& nu)
        :   IMEX<SweeperTrait, Enabled>()
          , _t0(0.0)
          , _nu(nu)
          , _lap(ndofs)
      {
        this->encap_factory()->set_size(ndofs * ndofs * ndofs);

        for (size_t yi = 0; yi < ndofs; ++yi) {
          this->_lap[yi] = vector<spatial_type>(ndofs);

          const spatial_type kyi = two_pi<spatial_type>()
                                   * ((yi <= ndofs / 2) ? spatial_type(yi) : spatial_type(yi) - spatial_type(ndofs));
          const spatial_type kyi_sqt = pfasst::almost_zero(pow(kyi, 2)) ? 0.0 : -pow(kyi, 2);

          for (size_t xi = 0; xi < ndofs; ++xi) {
            this->_lap[yi][xi] = vector<spatial_type>(ndofs);

            const spatial_type kxi = two_pi<spatial_type>()
                                     * ((xi <= ndofs / 2) ? spatial_type(xi) : spatial_type(xi) - spatial_type(ndofs));
            const spatial_type kxi_sqt = pfasst::almost_zero(pow(kxi, 2)) ? 0.0 : -pow(kxi, 2);

            for (size_t zi = 0; zi < ndofs; ++zi) {
              const spatial_type kzi = two_pi<spatial_type>()
                                       * ((zi <= ndofs / 2) ? spatial_type(zi) : spatial_type(zi) - spatial_type(ndofs));
              const spatial_type kzi_sqt = pfasst::almost_zero(pow(kzi, 2)) ? 0.0 : -pow(kzi, 2);

              this->_lap[yi][xi][zi] = kzi_sqt + kyi_sqt + kxi_sqt;
            }
          }
        }
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat3D<SweeperTrait, Enabled>::set_options()
      {
        IMEX<SweeperTrait, Enabled>::set_options();

        this->_nu = config::get_value<typename traits::spatial_type>("nu", 0.2);
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      Heat3D<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_type& t)
      {
        auto result = this->get_encap_factory().create();

        const size_t dofs_p_dim = cbrt(this->get_num_dofs());
        const spatial_type dx = 1.0 / spatial_type(dofs_p_dim);
        const spatial_type dy = 1.0 / spatial_type(dofs_p_dim);
        const spatial_type dz = 1.0 / spatial_type(dofs_p_dim);

        for (size_t yi = 0; yi < dofs_p_dim; ++yi) {
          for (size_t xi = 0; xi < dofs_p_dim; ++xi) {
            for (size_t zi = 0; zi < dofs_p_dim; ++zi) {
              result->data()[yi * dofs_p_dim + xi * dofs_p_dim + zi] = (sin(two_pi<spatial_type>() * yi * dy) + sin(two_pi<spatial_type>() * xi * dx) + sin(two_pi<spatial_type>() * zi * dz)) * exp(-t * 4 * pi_sqr<spatial_type>() * this->_nu);
            }
          }
        }

        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "EXACT t=" << t// << ": " << LOG_FLOAT << to_string(result)
        );

        return result;
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat3D<SweeperTrait, Enabled>::post_step()
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
      Heat3D<SweeperTrait, Enabled>::converged(const bool& pre_check)
      {
        const bool converged = IMEX<SweeperTrait, Enabled>::converged(pre_check);

        if (!pre_check) {
          assert(this->get_status() != nullptr);
          const time_type t = this->get_status()->get_time();
          const time_type dt = this->get_status()->get_dt();

  //         auto error = this->compute_error(t);
  //         auto rel_error = this->compute_relative_error(error, t);

          assert(this->get_quadrature() != nullptr);
          auto nodes = this->get_quadrature()->get_nodes();
          const auto num_nodes = this->get_quadrature()->get_num_nodes();
          nodes.insert(nodes.begin(), time_type(0.0));

          ML_CVLOG(1, this->get_logger_id(), "Observables after " << ((this->get_status()->get_iteration() == 0) ? string("prediction") : string("iteration ") + to_string(this->get_status()->get_iteration())));
          for (size_t m = 0; m < num_nodes; ++m) {
            ML_CVLOG(1, this->get_logger_id(), "  t["<<m<<"]=" << LOG_FIXED << (t + dt * nodes[m])
                                            << "      |abs residual| = " << LOG_FLOAT << this->_abs_res_norms[m]
                                            << "      |rel residual| = " << LOG_FLOAT << this->_rel_res_norms[m]
  //                                           << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[m])
  //                                           << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[m])
                    );
          }
          ML_CLOG(INFO, this->get_logger_id(), "  t["<<num_nodes<<"]=" << LOG_FIXED << (t + dt * nodes[num_nodes])
                                            << "      |abs residual| = " << LOG_FLOAT << this->_abs_res_norms[num_nodes]
                                            << "      |rel residual| = " << LOG_FLOAT << this->_rel_res_norms[num_nodes]
  //                                           << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[num_nodes])
  //                                           << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[num_nodes])
                 );
        }
        return converged;
      }

      template<class SweeperTrait, typename Enabled>
      size_t
      Heat3D<SweeperTrait, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory().size();
      }


      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_type>>
      Heat3D<SweeperTrait, Enabled>::compute_error(const typename SweeperTrait::time_type& t)
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
          ML_CVLOG(3, this->get_logger_id(), LOG_FIXED << "error t=" << (t + ds)// << ": "
//                                           << LOG_FLOAT << to_string(error[m])
                  );
        }

        return error;
      }

      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_type>>
      Heat3D<SweeperTrait, Enabled>::compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_type>>& error,
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
      Heat3D<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_type& t,
                                                       const shared_ptr<typename SweeperTrait::encap_type> u)
      {
        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "evaluating EXPLICIT part at t=" << t);
//         ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tu:   " << to_string(u));

        auto result = this->get_encap_factory().create();

        // taken form pySDC
        //   # xvalues = np.array([(i+1)*self.dx for i in range(self.nvars)])
        //   fexpl.values = np.zeros(self.nvars)  # -np.sin(np.pi * xvalues) * (np.sin(t) - self.nu * np.pi**2 * np.cos(t))
//         const spatial_type PI = pi<spatial_type>();
//         const spatial_type PIsqr = pi_sqr<spatial_type>();
//         const spatial_type dx = 1.0 / (spatial_type(this->get_num_dofs()) + 1);
//         for (size_t i = 0; i < this->get_num_dofs(); ++i) {
//           result->data()[i] = -1.0 * sin(PI * (i + 1) * dx) * (sin(t) - this->_nu * PIsqr * cos(t));
//         }
        result->zero();

        this->_num_expl_f_evals++;

//         ML_CVLOG(5, this->get_logger_id(), "\t  -> " << to_string(result));
        return result;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      Heat3D<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                       const shared_ptr<typename SweeperTrait::encap_type> u)
      {
        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "evaluating IMPLICIT part at t=" << t);
//         ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tu:   " << to_string(u));

        const spatial_type c = this->_nu / spatial_type(this->get_num_dofs());
        const size_t dofs_p_dim = cbrt(this->get_num_dofs());

        auto* z = this->_fft.forward(u);

        for (size_t yi = 0; yi < dofs_p_dim; ++yi) {
          for (size_t xi = 0; xi < dofs_p_dim; ++xi) {
            for (size_t zi = 0; zi < dofs_p_dim; ++zi) {
              z[yi * dofs_p_dim + xi * dofs_p_dim + zi] *= c * this->_lap[yi][xi][zi];
            }
          }
        }

        auto result = this->get_encap_factory().create();
        this->_fft.backward(result);

        this->_num_impl_f_evals++;

//         ML_CVLOG(5, this->get_logger_id(), "\t  -> " << to_string(result));
        return result;
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat3D<SweeperTrait, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_type> f,
                                                    shared_ptr<typename SweeperTrait::encap_type> u,
                                                    const typename SweeperTrait::time_type& t,
                                                    const typename SweeperTrait::time_type& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_type> rhs)
      {
        ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);
//         ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tf:   " << to_string(f));
//         ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\tu:   " << to_string(u));
//         ML_CVLOG(5, this->get_logger_id(), LOG_FLOAT << "\trhs: " << to_string(rhs));

        const spatial_type c = this->_nu * dt;
        const size_t dofs_p_dim = cbrt(this->get_num_dofs());

        auto* z = this->_fft.forward(rhs);

        for (size_t yi = 0; yi < dofs_p_dim; ++yi) {
          for (size_t xi = 0; xi < dofs_p_dim; ++xi) {
            for (size_t zi = 0; zi < dofs_p_dim; ++zi) {
              z[yi * dofs_p_dim + xi * dofs_p_dim + zi] /= (1.0 - c * this->_lap[yi][xi][zi]) * spatial_type(this->get_num_dofs());
            }
          }
        }

        this->_fft.backward(u);

        for (size_t yi = 0; yi < dofs_p_dim; ++yi) {
          for (size_t xi = 0; xi < dofs_p_dim; ++xi) {
            for (size_t zi = 0; zi < dofs_p_dim; ++zi) {
              f->data()[yi * dofs_p_dim + xi * dofs_p_dim + zi] = (u->get_data()[yi * dofs_p_dim + xi * dofs_p_dim + zi] - rhs->get_data()[yi * dofs_p_dim + xi * dofs_p_dim + zi]) / dt;
            }
          }
        }

        this->_num_impl_solves++;

//         ML_CVLOG(5, this->get_logger_id(), "\t->");
//         ML_CVLOG(5, this->get_logger_id(), "\t  f: " << to_string(f));
//         ML_CVLOG(5, this->get_logger_id(), "\t  u: " << to_string(u));
      }
    }  // ::pfasst::examples::heat1d
  }  // ::pfasst::examples
}  // ::pfasst
