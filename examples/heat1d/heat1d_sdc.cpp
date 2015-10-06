#include <memory>
using namespace std;

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>

#include <pfasst/encap/vector.hpp>
#include <pfasst/encap/eigen3_vector.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/contrib/spectral_transfer.hpp>

#include "heat1d_sweeper.hpp"

using encap_traits_t = pfasst::encap::vector_encap_traits<double, double, 1>;
// using encap_traits_t = pfasst::eigen3_encap_traits<double, double, 1>;


namespace pfasst
{
  namespace examples
  {
    namespace heat1d
    {
      using sweeper_t = Heat1D<pfasst::sweeper_traits<encap_traits_t>>;
      using pfasst::transfer_traits;
      using pfasst::contrib::SpectralTransfer;
      using pfasst::SDC;
      using pfasst::quadrature::QuadratureType;
      using heat1d_sdc_t = SDC<SpectralTransfer<transfer_traits<sweeper_t, sweeper_t, 1>>>;

      shared_ptr<heat1d_sdc_t> run_sdc(const size_t ndofs, const size_t nnodes,
                                               const QuadratureType& quad_type, const double& t_0,
                                               const double& dt, const double& t_end,
                                               const size_t niter)
      {
        using pfasst::quadrature::quadrature_factory;

        auto sdc = make_shared<heat1d_sdc_t>();

        auto sweeper = make_shared<sweeper_t>(ndofs);
        sweeper->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        sdc->add_sweeper(sweeper);
        sdc->set_options();

        sdc->status()->time() = t_0;
        sdc->status()->dt() = dt;
        sdc->status()->t_end() = t_end;
        sdc->status()->max_iterations() = niter;

        sdc->setup();

        sweeper->initial_state() = sweeper->exact(sdc->get_status()->get_time());

        sdc->run();
        sdc->post_run();

        return sdc;
      }
    }  // ::pfasst::examples::heat1d
  } // ::pfasst::examples
}  // ::pfasst


#ifndef PFASST_UNIT_TESTING
  int main(int argc, char** argv)
  {
    using pfasst::config::get_value;
    using pfasst::quadrature::QuadratureType;
    using pfasst::examples::heat1d::Heat1D;

    using sweeper_t = Heat1D<pfasst::sweeper_traits<encap_traits_t>>;

    pfasst::init(argc, argv, sweeper_t::init_opts);

    const size_t ndofs = get_value<size_t>("num_dofs", 8);
    const size_t nnodes = get_value<size_t>("num_nodes", 3);
    const QuadratureType quad_type = QuadratureType::GaussRadau;
    const double t_0 = 0.0;
    const double dt = get_value<double>("dt", 0.01);
    double t_end = get_value<double>("tend", -1);
    size_t nsteps = get_value<size_t>("num_steps", 0);
    if (t_end == -1 && nsteps == 0) {
      ML_CLOG(ERROR, "USER", "Either t_end or num_steps must be specified.");
      throw runtime_error("either t_end or num_steps must be specified");
    } else if (t_end != -1 && nsteps != 0) {
      if (!pfasst::almost_equal(t_0 + nsteps * dt, t_end)) {
        ML_CLOG(ERROR, "USER", "t_0 + nsteps * dt != t_end ("
                            << t_0 << " + " << nsteps << " * " << dt << " = " << (t_0 + nsteps * dt)
                            << " != " << t_end << ")");
        throw runtime_error("t_0 + nsteps * dt != t_end");
      }
    } else if (nsteps != 0) {
      t_end = t_0 + dt * nsteps;
    }
    const size_t niter = get_value<size_t>("num_iters", 5);

    pfasst::examples::heat1d::run_sdc(ndofs, nnodes, quad_type, t_0, dt, t_end, niter);
  }
#endif
