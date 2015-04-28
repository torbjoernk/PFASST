#include <memory>
using namespace std;

#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace ::testing;

#ifdef WITH_MPI
  #include <mpi.h>
#endif

#include <pfasst.hpp>
#ifdef WITH_MPI
  #include <pfasst/mpi_communicator.hpp>
#endif

#include "../examples/boris/boris_sweeper.hpp"
#include "../examples/boris/particle_cloud.hpp"
using namespace pfasst::examples::boris;

#ifdef WITH_MPI
TEST(CommunicationTest, BlockingSendRecv)
{
  pfasst::mpi::MPICommunicator comm(MPI_COMM_WORLD);
  ParticleCloudFactory<double> factory(2, 3, 1.0, 1.0);
  shared_ptr<Particle<double>> center = make_shared<Particle<double>>();
  center->pos()[0] = 1;
  center->vel()[0] = 2;
  shared_ptr<ParticleCloud<double>> test = dynamic_pointer_cast<ParticleCloud<double>>(factory.create(pfasst::encap::solution));
  shared_ptr<ParticleCloud<double>> copy = dynamic_pointer_cast<ParticleCloud<double>>(factory.create(pfasst::encap::solution));

  if (comm.rank() == comm.size() - 1) {
    test->distribute_around_center(center, comm);
    copy->copy(test);
    copy->broadcast(&comm);
  }
  if (comm.rank() % 2 == 0) {
    test->send(&comm, 0, true);
    test->recv(&comm, 0, true);
  } else if (comm.rank() % 2 == 1) {
    test->recv(&comm, 0, true);
    test->send(&comm, 0, true);
  }
  EXPECT_THAT(test->positions(), Pointwise(Eq(), copy->positions()));
  EXPECT_THAT(test->velocities(), Pointwise(Eq(), copy->velocities()));
}
#else
TEST(DISABLED_CommunicationTest, BlockingSendRecv)
{}
#endif

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
#ifdef WITH_MPI
  MPI_Init(&argc, &argv);
#endif
  pfasst::init(argc, argv,
               pfasst::examples::boris::init_opts<>,
               pfasst::examples::boris::init_logs<>);
  int result = 1, max_result = 0;  // GTest return value 1 (failure), 0 (success)
  // cppcheck-suppress redundantAssignment
  result = RUN_ALL_TESTS();
#ifdef WITH_MPI
  MPI_Allreduce(&result, &max_result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  return max_result;
}
