#include <memory>
using namespace std;

#include <mpi.h>

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/encap/vector.hpp>
#include <pfasst/comm/mpi_p2p.hpp>
#include <pfasst/controller/status.hpp>
#include <pfasst/controller/two_level_pfasst.hpp>
#include <pfasst/contrib/spectral_transfer.hpp>

#include "heat3d_space_parallel_sweeper.hpp"

using pfasst::encap::VectorEncapsulation;
using pfasst::quadrature::quadrature_factory;
using pfasst::quadrature::QuadratureType;
using pfasst::contrib::SpectralTransfer;
using pfasst::TwoLevelPfasst;
using temporal_comm_t = pfasst::comm::MpiP2P<pfasst::comm::temporal_communicator_tag>;
using spatial_comm_t = pfasst::comm::MpiP2P<pfasst::comm::spatial_communicator_tag>;

using pfasst::examples::heat3d::Heat3D;

using encap_t = VectorEncapsulation<double, double, 3>;
using sweeper_t = Heat3D<pfasst::space_parallel_sweeper_traits<typename encap_t::traits, spatial_comm_t>>;
using transfer_traits_t = pfasst::transfer_traits<sweeper_t, sweeper_t, 2>;
using spectral_t = SpectralTransfer<transfer_traits_t>;


namespace pfasst
{
  namespace examples
  {
    namespace heat3d
    {
      void run_pfasst(const size_t ndofs, const size_t nnodes, const QuadratureType& quad_type,
                      const double t_0, const double dt, const double t_end, const size_t niter,
                      const size_t np_space)
      {
        shared_ptr<temporal_comm_t> comm_time;
        shared_ptr<spatial_comm_t> comm_space;
        tie(comm_time, comm_space) = pfasst::comm::split_comm(np_space);

        TwoLevelPfasst<spectral_t, temporal_comm_t> pfasst;
        pfasst.communicator() = comm_time;

        auto coarse = make_shared<sweeper_t>(ndofs);
        coarse->spatial_communicator() = comm_space;
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = make_shared<sweeper_t>(ndofs * 2);
        fine->spatial_communicator() = comm_space;
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        auto transfer = make_shared<spectral_t>();

        pfasst.add_sweeper(coarse, true);
        pfasst.add_sweeper(fine, false);
        pfasst.add_transfer(transfer);
        pfasst.set_options();

        pfasst.status()->time() = t_0;
        pfasst.status()->dt() = dt;
        pfasst.status()->t_end() = t_end;
        pfasst.status()->max_iterations() = niter;

        pfasst.setup();

        coarse->initial_state() = coarse->exact(pfasst.get_status()->get_time());
        fine->initial_state() = fine->exact(pfasst.get_status()->get_time());

        pfasst.run();
        pfasst.post_run();
      }
    }  // ::pfasst::examples::heat2d
  } // ::pfasst::examples
}  // ::pfasst


int main(int argc, char** argv)
{
  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;

  MPI_Init(&argc, &argv);

  pfasst::init(argc, argv, sweeper_t::init_opts);
  pfasst::Status<double>::create_mpi_datatype();

  const size_t ndofs = get_value<size_t>("num_dofs", 4);
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

  const size_t np_space = get_value<size_t>("np_space", 1);

  pfasst::examples::heat3d::run_pfasst(ndofs, nnodes, quad_type, t_0, dt, t_end, niter, np_space);

  pfasst::Status<double>::free_mpi_datatype();

  MPI_Finalize();

  return 0;
}
