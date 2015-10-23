#include "fixtures/test_helpers.hpp"
using ::testing::Eq;
using ::testing::Lt;
using ::testing::Pointwise;

#include <vector>
using namespace std;

#include <pfasst/quadrature.hpp>
using pfasst::quadrature::QuadratureType;

#define PFASST_UNIT_TESTING
#include "examples/heat2d/heat2d_sdc.cpp"
#include "examples/heat2d/heat2d_mlsdc.cpp"


class SDCConvergence
  : public ::testing::Test
{};

TEST_F(SDCConvergence, single_step)
{
  ASSERT_FALSE(pfasst::config::has_value("abs_res_tol"));
  vector<double> abs_res(8);

  for (size_t k = 0; k < abs_res.size(); ++k) {
    auto sdc = pfasst::examples::heat2d::run_sdc(8, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.05, k);
    abs_res[k] = sdc->get_status()->get_abs_res_norm();
  }

  vector<double> expected{4.54e-2, 2.88e-3, 1.67e-4, 1.03e-5, 5.84e-7, 2.90e-8, 1.08e-9, 1.54e-10};
  EXPECT_THAT(abs_res, Pointwise(Lt(), expected));
}

TEST_F(SDCConvergence, two_step)
{
  ASSERT_FALSE(pfasst::config::has_value("abs_res_tol"));
  vector<double> abs_res(8);

  for (size_t k = 0; k < abs_res.size(); ++k) {
    auto sdc = pfasst::examples::heat2d::run_sdc(8, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.1, k);
    abs_res[k] = sdc->get_status()->get_abs_res_norm();
  }

  vector<double> expected{3.15e-2, 1.95e-3, 1.12e-4, 6.88e-6, 3.93e-7, 1.95e-8, 7.27e-10, 1.04e-10};
  EXPECT_THAT(abs_res, Pointwise(Lt(), expected));
}

TEST_F(SDCConvergence, abs_res_convergence)
{
  ConfigModder config;
  config.update<double>("abs_res_tol", 1e-10);
  ASSERT_TRUE(pfasst::config::has_value("abs_res_tol"));
  ASSERT_THAT(pfasst::config::get_value<double>("abs_res_tol", -1), Eq(1e-10));

  auto sdc = pfasst::examples::heat2d::run_sdc(8, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.05, 10);

  EXPECT_THAT(sdc->get_status()->get_iteration(), Eq(8));
}



class TwoLevelMLSDCConvergence
  : public ::testing::Test
{};

TEST_F(TwoLevelMLSDCConvergence, single_step)
{
  ASSERT_FALSE(pfasst::config::has_value("abs_res_tol"));
  vector<double> abs_res(6);

  for (size_t k = 0; k < abs_res.size(); ++k) {
    auto sdc = pfasst::examples::heat2d::run_mlsdc(8, 2, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.05, k);
    LOG(INFO) << to_string(sdc->get_status());
    abs_res[k] = sdc->get_status()->get_abs_res_norm();
  }

  vector<double> expected{2.88e-3, 1.03e-5, 2.90e-8, 1.54e-10, 1.30e-12, 7.10e-15};
  EXPECT_THAT(abs_res, Pointwise(Lt(), expected));
}

TEST_F(TwoLevelMLSDCConvergence, two_step)
{
  ASSERT_FALSE(pfasst::config::has_value("abs_res_tol"));
  vector<double> abs_res(6);

  for (size_t k = 0; k < abs_res.size(); ++k) {
    auto sdc = pfasst::examples::heat2d::run_mlsdc(8, 2, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.1, k);
    abs_res[k] = sdc->get_status()->get_abs_res_norm();
  }

  vector<double> expected{1.95e-3, 6.88e-6, 1.95e-8, 1.04e-10, 8.71e-13, 4.85e-15};
  EXPECT_THAT(abs_res, Pointwise(Lt(), expected));
}

TEST_F(TwoLevelMLSDCConvergence, abs_res_convergence)
{
  ConfigModder config;
  config.update<double>("abs_res_tol", 1e-10);
  ASSERT_THAT(pfasst::config::get_value<double>("abs_res_tol", -1), Eq(1e-10));

  auto sdc = pfasst::examples::heat2d::run_mlsdc(8, 2, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.05, 10);

  EXPECT_THAT(sdc->get_status()->get_iteration(), Eq(4));
}


TEST_MAIN()
