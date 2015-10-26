#include <memory>
#include <stdexcept>
using std::shared_ptr;

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/encap/vector.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/contrib/spectral_transfer.hpp>

#include "heat3d_sweeper.hpp"

using pfasst::quadrature::QuadratureType;

using DataVectorTraits = pfasst::encap::vector_encap_traits<double, double, 3>;


namespace pfasst
{
  namespace examples
  {
    namespace heat3d
    {
      shared_ptr<
        pfasst::SDC<
          pfasst::contrib::SpectralTransfer<
            pfasst::transfer_traits<
              Heat3D<pfasst::sweeper_traits<DataVectorTraits>>,
              Heat3D<pfasst::sweeper_traits<DataVectorTraits>>,
              1
            >
          >
        >
      >
      run_sdc(const size_t& ndofs, const size_t& nnodes, const QuadratureType& quad_type,
              const double& t_0, const double& dt, const double& t_end, const size_t& niter)
      {
        using pfasst::quadrature::quadrature_factory;
        using pfasst::contrib::SpectralTransfer;
        using pfasst::SDC;

        using SweeperType = Heat3D<pfasst::sweeper_traits<DataVectorTraits>>;
        using TransferType = SpectralTransfer<pfasst::transfer_traits<SweeperType, SweeperType, 1>>;

        auto sdc = std::make_shared<SDC<TransferType>>();

        auto sweeper = std::make_shared<SweeperType>(ndofs);
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
    }  // ::pfasst::examples::heat2d
  } // ::pfasst::examples
}  // ::pfasst


#ifndef PFASST_UNIT_TESTING
  int main(int argc, char** argv)
  {
    using pfasst::config::get_value;
    using pfasst::quadrature::QuadratureType;
    using pfasst::examples::heat3d::Heat3D;

    using SweeperType = Heat3D<pfasst::sweeper_traits<DataVectorTraits>>;

    pfasst::init(argc, argv, SweeperType::init_opts);

    const size_t ndofs = get_value<size_t>("num_dofs", 2);
    const size_t nnodes = get_value<size_t>("num_nodes", 3);
    const QuadratureType quad_type = QuadratureType::GaussRadau;
    const double t_0 = 0.0;
    const double dt = get_value<double>("dt", 0.01);
    double t_end = get_value<double>("tend", -1);
    size_t nsteps = get_value<size_t>("num_steps", 0);
    if (t_end == -1 && nsteps == 0) {
      ML_CLOG(ERROR, "USER", "Either t_end or num_steps must be specified.");
      throw std::runtime_error("either t_end or num_steps must be specified");
    } else if (t_end != -1 && nsteps != 0) {
      if (!pfasst::almost_equal(t_0 + nsteps * dt, t_end)) {
        ML_CLOG(ERROR, "USER", "t_0 + nsteps * dt != t_end ("
                            << t_0 << " + " << nsteps << " * " << dt << " = " << (t_0 + nsteps * dt)
                            << " != " << t_end << ")");
        throw std::runtime_error("t_0 + nsteps * dt != t_end");
      }
    } else if (nsteps != 0) {
      t_end = t_0 + dt * nsteps;
    }
    const size_t niter = get_value<size_t>("num_iters", 5);

    pfasst::examples::heat3d::run_sdc(ndofs, nnodes, quad_type, t_0, dt, t_end, niter);
  }
#endif
