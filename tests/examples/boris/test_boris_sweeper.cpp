#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace ::testing;

#include "pfasst/logging.hpp"

#define PFASST_UNIT_TESTING
#include "../examples/boris/boris_sweeper.hpp"
#include "../examples/boris/boris_sdc.cpp"
#undef PFASST_UNIT_TESTING
using namespace pfasst::examples::boris;


TEST(EnergyDriftAndResidual, SingleStep)
{
  const size_t num_iter = 9;
  auto errors_map = run_boris_sdc<double>(1, 0.015625, 5, 1, num_iter+1, 0.0, 0.0);
  ASSERT_THAT(errors_map, SizeIs(num_iter+1));

  auto final_error = errors_map.rbegin()->second;

  EXPECT_THAT(final_error.e_drift, DoubleNear(0.0, 2e-12));
  EXPECT_THAT(final_error.res, DoubleNear(0.0, 1.5e-14));
}

TEST(EnergyDriftAndResidual, MultiStep)
{
  const size_t num_iter = 9;
  const size_t num_steps = 10;
  auto errors_map = run_boris_sdc<double>(num_steps, 0.015625, 5, 1, num_iter+1, 0.0, 0.0);
  ASSERT_THAT(errors_map, SizeIs((num_iter+1) * num_steps));

  auto final_error = errors_map.rbegin()->second;

  EXPECT_THAT(final_error.e_drift, DoubleNear(0.0, 1.1e-11));
  EXPECT_THAT(final_error.res, DoubleNear(0.0, 1.5e-14));
}


int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
#ifdef WITH_MPI
  MPI_Init(&argc, &argv);
#endif
  pfasst::init(argc, argv,
               pfasst::examples::boris::init_opts<>,
               pfasst::examples::boris::init_logs<>);
  int result = 1, max_result;  // GTest return value 1 (failure), 0 (success)
  result = RUN_ALL_TESTS();
#ifdef WITH_MPI
  MPI_Allreduce(&result, &max_result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return max_result;
#else
  return result;
#endif
}
