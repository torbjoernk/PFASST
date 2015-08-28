#include "fixtures/test_helpers.hpp"

#include <vector>
using namespace std;

#include <pfasst/quadrature.hpp>
#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
using pfasst::quadrature::QuadratureType;
typedef pfasst::vector_encap_traits<double, double>                           EncapTraits;

#include "examples/heat1d/heat1d_sweeper.hpp"
typedef pfasst::examples::heat1d::Heat1D<pfasst::sweeper_traits<EncapTraits>> SweeperType;

#define PFASST_UNIT_TESTING
#include "examples/heat1d/heat1d_sdc.cpp"


class ProblemSetup
  : public ::testing::Test
{
  protected:
    typedef          SweeperType              sweeper_type;
    typedef typename sweeper_type::encap_type encap_type;

    shared_ptr<sweeper_type> sweeper;

    shared_ptr<pfasst::Status<double>> status = make_shared<pfasst::Status<double>>();

    vector<double> exact_t0{   0.000000000000000000000000
                            ,  0.7071067811865474617150085
                            ,  1.000000000000000000000000
                            ,  0.7071067811865475727373109
                            ,  1.224646799147353207173764e-16
                            , -0.7071067811865474617150085
                            , -1.000000000000000000000000
                            , -0.7071067811865476837596134};
    vector<double> exact_t001{   0.000000000000000000000000e+00
                              ,  0.7015456730923740336081096
                              ,  0.9921354055113971170953846
                              ,  0.7015456730923741446304120
                              ,  1.215015448680293835571449e-16
                              , -0.7015456730923740336081096
                              , -0.9921354055113971170953846
                              , -0.7015456730923742556527145};

    virtual void SetUp()
    {
      sweeper = make_shared<sweeper_type>(8);
      sweeper->status() = status;
    }
};

TEST_F(ProblemSetup, computes_exact_solution_at_t0)
{
  auto exact = sweeper->exact(0.0);

  EXPECT_THAT(exact->get_data(), Pointwise(DoubleNear(), exact_t0));
}

TEST_F(ProblemSetup, computes_exact_solution_at_t01)
{
  auto exact = sweeper->exact(0.01);

  EXPECT_THAT(exact->get_data(), Pointwise(DoubleNear(), exact_t001));
}


class Convergence
  : public ::testing::Test
{};

TEST_F(Convergence, single_step)
{
  vector<double> abs_res(8);

  for (size_t k = 0; k < abs_res.size(); ++k) {
    auto sdc = pfasst::examples::heat1d::run_sdc(8, 3, QuadratureType::GaussRadau, 0.0, 0.05, 0.05, k+1);
    abs_res[k] = sdc->get_status()->get_abs_res_norm();
  }

  vector<double> expected{1.44e-3, 8.31e-5, 5.11e-6, 2.92e-7, 1.45e-8, 5.39e-10, 7.70e-11, 7.66e-12};
  EXPECT_THAT(abs_res, Pointwise(Lt(), expected));
}


TEST_MAIN()
