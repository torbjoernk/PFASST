#include "fixtures/test_helpers.hpp"
using ::testing::Eq;
using ::testing::NotNull;

#include <memory>
#include <stdexcept>
using std::shared_ptr;

#include <pfasst/controller/controller.hpp>
using pfasst::Controller;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>

#include <pfasst/transfer/traits.hpp>
#include <pfasst/transfer/polynomial.hpp>

#include "controller/mocks.hpp"
#include "sweeper/mocks.hpp"
#include "transfer/mocks.hpp"

using encap_traits_t = pfasst::encap::vector_encap_traits<double, double, 1>;
using encap_t = pfasst::encap::Encapsulation<encap_traits_t>;
using sweeper_t = SweeperMock<pfasst::sweeper_traits<encap_traits_t>>;
using transfer_traits_t = pfasst::transfer_traits<sweeper_t, sweeper_t, 2>;
using transfer_t = TransferMock<transfer_traits_t>;


using ControllerTypes = ::testing::Types<Controller<transfer_t>>;
INSTANTIATE_TYPED_TEST_CASE_P(Controller, Concepts, ControllerTypes);


class Interface
  : public ::testing::Test
{
  protected:
    shared_ptr<Controller<transfer_t>> controller;
    shared_ptr<pfasst::Status<double>> status;

    virtual void SetUp()
    {
      this->controller = std::make_shared<Controller<transfer_t>>();
      this->status = std::make_shared<pfasst::Status<double>>();
    }
};

TEST_F(Interface, has_a_status)
{
  EXPECT_THAT(controller->get_status(), NotNull());
}

TEST_F(Interface, status_can_be_assigned)
{
  controller->status() = status;
  EXPECT_THAT(controller->get_status(), Eq(status));
}

TEST_F(Interface, status_can_be_modified)
{
  controller->status()->time() = 42.0;
  EXPECT_THAT(controller->get_status()->get_time(), Eq(42.0));
}


class Setup
  : public ::testing::Test
{
  protected:
    shared_ptr<Controller<transfer_t>> controller;
    shared_ptr<pfasst::Status<double>> status;
    shared_ptr<transfer_t> transfer;

    virtual void SetUp()
    {
      this->controller = std::make_shared<Controller<transfer_t>>();
      this->status = std::make_shared<pfasst::Status<double>>();
    }
};

TEST_F(Setup, setup_required_for_running)
{
  controller->status()->t_end() = 4.2;
  controller->status()->dt() = 0.1;
  controller->status()->max_iterations() = 1;

  ASSERT_FALSE(controller->is_ready());
  EXPECT_THROW(controller->run(), std::logic_error);

  controller->setup();
  EXPECT_TRUE(controller->is_ready());

  controller->run();
}


class Logic
  : public ::testing::Test
{
  protected:
    shared_ptr<Controller<transfer_t>> controller;
    shared_ptr<pfasst::Status<double>> status;

    virtual void SetUp()
    {
      this->controller = std::make_shared<Controller<transfer_t>>();
      this->status = std::make_shared<pfasst::Status<double>>();
    }
};

TEST_F(Logic, advance_in_time_with_sufficient_t_end)
{
  controller->status()->dt() = 0.1;
  controller->status()->time() = 1.0;
  controller->status()->step() = 1;
  controller->status()->t_end() = 2.0;

  EXPECT_TRUE(controller->advance_time());
  EXPECT_THAT(controller->get_status()->get_time(), Eq(1.1));
  EXPECT_THAT(controller->get_status()->get_step(), Eq(2));
}

TEST_F(Logic, advance_in_time_with_insufficient_t_end)
{
  controller->status()->dt() = 0.1;
  controller->status()->time() = 1.0;
  controller->status()->step() = 1;
  controller->status()->t_end() = 1.0;

  EXPECT_FALSE(controller->advance_time());
  EXPECT_THAT(controller->get_status()->get_time(), Eq(1.0));
  EXPECT_THAT(controller->get_status()->get_step(), Eq(1));
}

TEST_F(Logic, advance_in_time_multiple_steps_at_once)
{
  controller->status()->dt() = 0.1;
  controller->status()->time() = 1.0;
  controller->status()->step() = 1;
  controller->status()->t_end() = 2.0;

  EXPECT_TRUE(controller->advance_time(3));
  EXPECT_THAT(controller->get_status()->get_time(), Eq(1.3));
  EXPECT_THAT(controller->get_status()->get_step(), Eq(4));
}


TEST_F(Logic, advance_iteration_with_exceeding_max_iteration_threshold)
{
  controller->status()->iteration() = 1;
  controller->status()->max_iterations() = 1;
  ASSERT_THAT(controller->get_status()->get_iteration(), Eq(1));
  ASSERT_THAT(controller->get_status()->get_max_iterations(), Eq(1));

  EXPECT_FALSE(controller->advance_iteration());
  EXPECT_THAT(controller->get_status()->get_iteration(), Eq(1));
}

TEST_F(Logic, advance_iteration)
{
  controller->status()->iteration() = 1;
  controller->status()->max_iterations() = 5;
  ASSERT_THAT(controller->get_status()->get_iteration(), Eq(1));
  ASSERT_THAT(controller->get_status()->get_max_iterations(), Eq(5));

  EXPECT_TRUE(controller->advance_iteration());
  EXPECT_THAT(controller->get_status()->get_iteration(), Eq(2));
}


TEST_MAIN()
