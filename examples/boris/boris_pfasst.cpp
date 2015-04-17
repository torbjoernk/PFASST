#include <memory>

#include <pfasst.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/controller/pfasst.hpp>
#include <pfasst/mpi_communicator.hpp>
using namespace pfasst;

#include "particle.hpp"
#include "particle_cloud.hpp"
#include "bindings/wrapper_interface.hpp"
#include "bindings/wrapper_simple_physics_solver.hpp"
#include "boris_sweeper.hpp"
#include "injective_transfer.hpp"


namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      template<typename scalar>
      error_map<scalar> run_boris_pfasst(const size_t nsteps, const scalar dt, const size_t nnodes,
                                         const size_t nparticles, const size_t niters,
                                         const double abs_res_tol, const double rel_res_tol,
                                         const size_t nblocks)
      {
        int size_world = -1, rank_world = -1;
        MPI_Comm_size(MPI_COMM_WORLD, &size_world);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
        MPI_Comm BORIS_COMM_SPACE;
        MPI_Comm BORIS_COMM_TIME;

        CLOG(DEBUG, "default") << "splitting up communicators into " << nblocks << " time blocks";

        int color_time = rank_world / (size_world / nblocks);
        int err_split_time = MPI_Comm_split(MPI_COMM_WORLD, color_time, rank_world, &BORIS_COMM_TIME);
        if (err_split_time != MPI_SUCCESS) {
          throw pfasst::mpi::MPIError("failed to split MPI_COMM_WORLD into time blocks");
        }
        int size_time = -1, rank_time = -1;
        MPI_Comm_size(BORIS_COMM_TIME, &size_time);
        MPI_Comm_rank(BORIS_COMM_TIME, &rank_time);
        CLOG(DEBUG, "default") << "got time rank " << rank_time << " of " << size_time << " ranks in time block " << color_time;

        int color_space = rank_world % size_time;
        int err_split_space = MPI_Comm_split(MPI_COMM_WORLD, color_space, rank_world, &BORIS_COMM_SPACE);
        if (err_split_space != MPI_SUCCESS) {
          throw pfasst::mpi::MPIError("failed to split MPI_COMM_WORLD into space blocks");
        }
        int size_space = -1, rank_space = -1;
        MPI_Comm_size(BORIS_COMM_SPACE, &size_space);
        MPI_Comm_rank(BORIS_COMM_SPACE, &rank_space);
        CLOG(DEBUG, "default") << "got space rank " << rank_space << " of " << size_space << " ranks in space block " << color_space;

        PFASST<> controller;
        mpi::MPICommunicator comm_time(BORIS_COMM_TIME);
        mpi::MPICommunicator comm_space(BORIS_COMM_SPACE);
        controller.set_comm(&comm_time);

        const double mass = 1.0;
        const double charge = 1.0;

        shared_ptr<bindings::WrapperInterface<double, double>> impl_solver = \
          make_shared<bindings::WrapperSimplePhysicsSolver<double, double>>(BORIS_COMM_SPACE);
        bindings::setup(dynamic_pointer_cast<bindings::WrapperSimplePhysicsSolver<double, double>>(impl_solver));

        if (nparticles % size_space != 0) {
          CLOG(ERROR, "Boris") << "Total number of particles (" << nparticles << ") "
                               << "must be a multiple of number of ranks in space blocks (" << size_space << ")";
          throw pfasst::ValueError("number of particles must be multiple of ranks in space blocks");
        }

        // fine level
        auto quad1        = quadrature::quadrature_factory<double>(nnodes,
                                                                   quadrature::QuadratureType::GaussLobatto);
        auto factory1     = make_shared<ParticleCloudFactory<double>>(nparticles / size_space, 3, mass, charge);
        string data_file1 = "s" + to_string(nsteps) + "_i" + to_string(niters)
                            + "_dt" + to_string(dt) + "_m" + to_string(nnodes)
                            + "_p" + to_string(nparticles)
                            + "_np" + to_string(size_world) + "-" + to_string(rank_world)
                            + "_level1.csv";
        auto sweeper1     = make_shared<BorisSweeper<double, double>>(impl_solver, data_file1);
        auto transfer1    = make_shared<InjectiveTransfer<double, double>>();
        sweeper1->set_quadrature(quad1);
        sweeper1->set_factory(factory1);
        sweeper1->set_residual_tolerances(abs_res_tol, rel_res_tol);
        controller.add_level(sweeper1, transfer1);

        // coarse level
        auto quad2        = quadrature::quadrature_factory<double>(nnodes,
                                                                   quadrature::QuadratureType::GaussLobatto);
        auto factory2     = make_shared<ParticleCloudFactory<double>>(nparticles / size_space, 3, mass, charge);
        string data_file2 = "s" + to_string(nsteps) + "_i" + to_string(niters)
                            + "_dt" + to_string(dt) + "_m" + to_string(nnodes)
                            + "_p" + to_string(nparticles)
                            + "_np" + to_string(size_world) + "-" + to_string(rank_world)
                            + "_level2.csv";
        auto sweeper2     = make_shared<BorisSweeper<double, double>>(impl_solver, data_file2);
        auto transfer2    = make_shared<InjectiveTransfer<double, double>>();
        sweeper2->set_quadrature(quad2);
        sweeper2->set_factory(factory2);
        controller.add_level(sweeper2, transfer2);

        controller.set_duration(0.0, nsteps*dt, dt, niters);
        controller.set_options();
        controller.setup();

        shared_ptr<Particle<double>> center = make_shared<Particle<double>>();
        center->pos()[0] = 10;
        center->vel()[0] = 100;
        center->vel()[2] = 100;

        auto fine_sweeper = controller.get_finest<BorisSweeper<double, double>>();
        shared_ptr<ParticleCloud<double>> q0 = \
          dynamic_pointer_cast<ParticleCloud<double>>(fine_sweeper->get_start_state());
        q0->distribute_around_center(center, comm_space);
        CLOG(INFO, "Boris") << "Initial Particles (fine) : "
                            << *(dynamic_pointer_cast<ParticleCloud<double>>(fine_sweeper->get_start_state()));
        fine_sweeper->set_initial_energy();

        controller.run();

        return fine_sweeper->get_errors();
      }
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  pfasst::init(argc, argv,
               pfasst::examples::boris::init_opts<>,
               pfasst::examples::boris::init_logs<>);

  const size_t nsteps     = pfasst::config::get_value<size_t>("num_steps", 1);
  const double dt         = pfasst::config::get_value<double>("delta_step", 0.015625);
  const size_t nnodes     = pfasst::config::get_value<size_t>("num_nodes", 5);
  const size_t nparticles = pfasst::config::get_value<size_t>("num_particles", 1);
  const size_t niters     = pfasst::config::get_value<size_t>("num_iter", 2);
  const double abs_res_tol = pfasst::config::get_value<double>("abs_res_tol", 0.0);
  const double rel_res_tol = pfasst::config::get_value<double>("rel_res_tol", 0.0);
  const size_t num_blocks = pfasst::config::get_value<size_t>("num_blocks", 1);

  CLOG(INFO, "Boris") << "nsteps=" << nsteps << ", "
                      << "dt=" << dt << ", "
                      << "nnodes=" << nnodes << ", "
                      << "nparticles=" << nparticles << ", "
                      << "niter=" << niters << ", "
                      << "abs res=" << abs_res_tol << ", "
                      << "rel res=" << rel_res_tol << ", "
                      << "nblocks=" << num_blocks;

  pfasst::examples::boris::run_boris_pfasst<double>(nsteps, dt, nnodes, nparticles, niters, abs_res_tol, rel_res_tol, num_blocks);
  MPI_Finalize();
}
#endif
