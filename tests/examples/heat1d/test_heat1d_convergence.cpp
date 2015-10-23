#include "fixtures/test_helpers.hpp"
using ::testing::Eq;
using ::testing::Lt;
using ::testing::Pointwise;

#include <vector>
using namespace std;

#include <pfasst/quadrature.hpp>
using pfasst::quadrature::QuadratureType;

#define PFASST_UNIT_TESTING
#include "examples/heat1d/heat1d_sdc.cpp"
#include "examples/heat1d/heat1d_mlsdc.cpp"


class SDCConvergence
  : public ::testing::Test
{};

TEST_F(SDCConvergence, single_step)
{
  ASSERT_FALSE(pfasst::config::has_value("abs_res_tol"));
  vector<double> abs_res(8);

  for (size_t k = 0; k < abs_res.size(); ++k) {
    auto sdc = pfasst::examples::heat1d::run_sdc(8, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.05, k+1);
    abs_res[k] = sdc->get_status()->get_abs_res_norm();
  }

  vector<double> expected{1.44e-3, 8.31e-5, 5.11e-6, 2.92e-7, 1.45e-8, 5.39e-10, 7.70e-11, 7.66e-12};
  EXPECT_THAT(abs_res, Pointwise(Lt(), expected));
}

TEST_F(SDCConvergence, two_step)
{
  ASSERT_FALSE(pfasst::config::has_value("abs_res_tol"));
  vector<double> abs_res(8);

  for (size_t k = 0; k < abs_res.size(); ++k) {
    auto sdc = pfasst::examples::heat1d::run_sdc(8, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.1, k+1);
    abs_res[k] = sdc->get_status()->get_abs_res_norm();
  }

  vector<double> expected{9.72e-4, 5.6e-5, 3.44e-6, 1.97e-7, 9.75e-9, 3.64e-10, 5.19e-11, 5.17e-12};
  EXPECT_THAT(abs_res, Pointwise(Lt(), expected));
}

TEST_F(SDCConvergence, abs_res_convergence)
{
  ConfigModder config;
  config.update<double>("abs_res_tol", 1e-10);
  ASSERT_TRUE(pfasst::config::has_value("abs_res_tol"));
  ASSERT_THAT(pfasst::config::get_value<double>("abs_res_tol", -1), Eq(1e-10));

  auto sdc = pfasst::examples::heat1d::run_sdc(8, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.05, 10);

  EXPECT_THAT(sdc->get_status()->get_iteration(), Eq(7));
}



class TwoLevelMLSDCConvergence
  : public ::testing::Test
{};

TEST_F(TwoLevelMLSDCConvergence, single_step)
{
  ASSERT_FALSE(pfasst::config::has_value("abs_res_tol"));
  vector<double> abs_res(6);

  for (size_t k = 0; k < abs_res.size(); ++k) {
    auto sdc = pfasst::examples::heat1d::run_mlsdc(8, 2, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.05, k);
    LOG(INFO) << to_string(sdc->get_status());
    abs_res[k] = sdc->get_status()->get_abs_res_norm();
  }

  vector<double> expected{1.44e-3, 5.11e-6, 1.45e-8, 7.70e-11, 6.46e-13, 3.70e-15};
  EXPECT_THAT(abs_res, Pointwise(Lt(), expected));
}

TEST_F(TwoLevelMLSDCConvergence, two_step)
{
  ASSERT_FALSE(pfasst::config::has_value("abs_res_tol"));
  vector<double> abs_res(6);

  for (size_t k = 0; k < abs_res.size(); ++k) {
    auto sdc = pfasst::examples::heat1d::run_mlsdc(8, 2, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.1, k);
    abs_res[k] = sdc->get_status()->get_abs_res_norm();
  }

  vector<double> expected{9.72e-4, 3.44e-6, 9.75e-9, 5.19e-11, 4.36e-13, 2.39e-15};
  EXPECT_THAT(abs_res, Pointwise(Lt(), expected));
}

TEST_F(TwoLevelMLSDCConvergence, abs_res_convergence)
{
  ConfigModder config;
  config.update<double>("abs_res_tol", 1e-10);
  ASSERT_THAT(pfasst::config::get_value<double>("abs_res_tol", -1), Eq(1e-10));

  auto sdc = pfasst::examples::heat1d::run_mlsdc(8, 2, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.05, 10);

  EXPECT_THAT(sdc->get_status()->get_iteration(), Eq(3));
}


TEST_MAIN()
