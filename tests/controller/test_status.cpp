#include "fixtures/test_helpers.hpp"
using ::testing::Eq;

#include <memory>
#include <stdexcept>
using std::make_shared;

#include <pfasst/controller/status.hpp>

#include "comm/mocks.hpp"
using comm_t = CommMock;

using StatusTypes = ::testing::Types<pfasst::Status<double>>;
INSTANTIATE_TYPED_TEST_CASE_P(Status, Concepts, StatusTypes);


class Interface
  : public ::testing::Test
{
  protected:
    pfasst::Status<double> status;
};

TEST_F(Interface, has_a_step)
{
  ASSERT_THAT(status.get_step(), Eq(0));

  status.step() = 1;
  EXPECT_THAT(status.get_step(), Eq(1));
}

TEST_F(Interface, has_an_iteration)
{
  ASSERT_THAT(status.get_iteration(), Eq(0));

  status.iteration() = 1;
  EXPECT_THAT(status.get_iteration(), Eq(1));
}

TEST_F(Interface, has_a_time_point)
{
  ASSERT_THAT(status.get_time(), Eq(0.0));

  status.time() = 1.42;
  EXPECT_THAT(status.get_time(), Eq(1.42));
}

TEST_F(Interface, has_a_time_delta)
{
  ASSERT_THAT(status.get_dt(), Eq(0.0));

  status.dt() = 0.42;
  EXPECT_THAT(status.get_dt(), Eq(0.42));
}

TEST_F(Interface, has_a_primary_state)
{
  ASSERT_THAT(status.get_primary_state(), Eq((+pfasst::PrimaryState::UNKNOWN_PRIMARY)));

  status.set_primary_state(pfasst::PrimaryState::CONVERGED);
  EXPECT_THAT(status.get_primary_state(), Eq((+pfasst::PrimaryState::CONVERGED)));
}

TEST_F(Interface, has_a_secondary_state)
{
  ASSERT_THAT(status.get_secondary_state(), Eq((+pfasst::SecondaryState::UNKNOWN_SECONDARY)));

  status.set_primary_state(pfasst::PrimaryState::ITERATING);
  status.set_secondary_state(pfasst::SecondaryState::ITER_FINE);
  EXPECT_THAT(status.get_secondary_state(), Eq((+pfasst::SecondaryState::ITER_FINE)));
}

TEST_F(Interface, combination_of_primary_and_secondary_states_are_validated)
{
  ASSERT_THAT(status.get_primary_state(), Eq((+pfasst::PrimaryState::UNKNOWN_PRIMARY)));
  ASSERT_THAT(status.get_secondary_state(), Eq((+pfasst::SecondaryState::UNKNOWN_SECONDARY)));

  status.set_primary_state(pfasst::PrimaryState::ITERATING);
  EXPECT_THROW(status.set_secondary_state(pfasst::SecondaryState::CONV_CHECK), std::runtime_error);
}

TEST_F(Interface, has_an_absolute_residual_norm)
{
  ASSERT_THAT(status.get_abs_res_norm(), Eq(0.0));

  status.abs_res_norm() = 0.1;
  EXPECT_THAT(status.get_abs_res_norm(), Eq(0.1));
}

TEST_F(Interface, has_a_relative_residual_norm)
{
  ASSERT_THAT(status.get_rel_res_norm(), Eq(0.0));

  status.rel_res_norm() = 0.1;
  EXPECT_THAT(status.get_rel_res_norm(), Eq(0.1));
}


class Communication
  : public ::testing::Test
{
  protected:
    shared_ptr<pfasst::Status<double>> status;
    shared_ptr<comm_t>                 comm;

  public:
    virtual void SetUp()
    {
      this->status = make_shared<pfasst::Status<double>>();
      this->comm = make_shared<comm_t>();
    }
};

TEST_F(Communication, can_be_send) {
  status->send(comm, 1, 0, true);

  status->send(comm, 1, 0, false);
}

TEST_F(Communication, can_be_received) {
  status->recv(comm, 1, 0, true);

  status->recv(comm, 1, 0, false);
}


TEST_MAIN()
