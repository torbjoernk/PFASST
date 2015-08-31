#include "fixtures/test_helpers.hpp"

#include <vector>
using namespace std;

#include <pfasst/quadrature.hpp>
using pfasst::quadrature::QuadratureType;

#define PFASST_UNIT_TESTING
#include "examples/heat1d/heat1d_sdc.cpp"


class SDCConvergence
  : public ::testing::Test
{};

TEST_F(SDCConvergence, single_step)
{
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
  vector<double> abs_res(8);
  
  for (size_t k = 0; k < abs_res.size(); ++k) {
    auto sdc = pfasst::examples::heat1d::run_sdc(8, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.1, k+1);
    abs_res[k] = sdc->get_status()->get_abs_res_norm();
  }
  
  vector<double> expected{9.72e-4, 5.6e-5, 3.44e-6, 1.97e-7, 9.75e-9, 3.64e-10, 5.19e-11, 5.17e-12};
  EXPECT_THAT(abs_res, Pointwise(Lt(), expected));
}


TEST_MAIN()
