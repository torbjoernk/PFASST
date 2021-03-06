#include "fixtures/test_helpers.hpp"
using ::testing::_;
using ::testing::Eq;
using ::testing::Mock;
using ::testing::NiceMock;
using ::testing::NotNull;
using ::testing::Return;
using ::testing::ReturnRef;

#include <memory>
#include <stdexcept>
#include <vector>
using std::shared_ptr;
using std::vector;

#include <pfasst/controller/sdc.hpp>
using pfasst::SDC;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>

#include <pfasst/transfer/traits.hpp>
#include <pfasst/transfer/polynomial.hpp>

#include "controller/mocks.hpp"
#include "quadrature/mocks.hpp"
#include "sweeper/mocks.hpp"
#include "transfer/mocks.hpp"

using encap_traits_t = pfasst::encap::vector_encap_traits<double, double, 1>;
using encap_t = pfasst::encap::Encapsulation<encap_traits_t>;
using sweeper_t = NiceMock<SweeperMock<pfasst::sweeper_traits<encap_traits_t>>>;
using transfer_traits_t = pfasst::transfer_traits<sweeper_t, sweeper_t, 1>;
using transfer_t = NiceMock<TransferMock<transfer_traits_t>>;
using quadrature_t = NiceMock<QuadratureMock<double>>;


using SDCTypes = ::testing::Types<SDC<transfer_t>>;
INSTANTIATE_TYPED_TEST_CASE_P(SDC, Concepts, SDCTypes);


class Interface
  : public ::testing::Test
{
  protected:
    shared_ptr<SDC<transfer_t>> controller;
    shared_ptr<pfasst::Status<double>> status;

    virtual void SetUp()
    {
      this->controller = std::make_shared<SDC<transfer_t>>();
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
  controller->status() = status;
  controller->status()->time() = 42.0;
  EXPECT_THAT(controller->get_status()->get_time(), Eq(42.0));
}


class Setup
  : public ::testing::Test
{
  protected:
    shared_ptr<SDC<transfer_t>> controller;

    vector<double>                     nodes{0.0, 0.5, 1.0};
    shared_ptr<pfasst::Status<double>> status;
    shared_ptr<sweeper_t>              sweeper;
    shared_ptr<transfer_t>             transfer;
    shared_ptr<quadrature_t>           quad;
    shared_ptr<encap_t>                encap;

    virtual void SetUp()
    {
      this->encap = std::make_shared<encap_t>(vector<double>{0.0, 0.0});
      this->controller = std::make_shared<SDC<transfer_t>>();
      this->status = std::make_shared<pfasst::Status<double>>();
      this->controller->status() = status;

      this->quad = std::make_shared<quadrature_t>();
      ON_CALL(*(this->quad.get()), right_is_node()).WillByDefault(Return(true));
      ON_CALL(*(this->quad.get()), get_nodes()).WillByDefault(ReturnRef(this->nodes));
      ON_CALL(*(this->quad.get()), get_num_nodes()).WillByDefault(Return(3));

      this->sweeper = std::make_shared<sweeper_t>();
      ON_CALL(*(this->sweeper.get()), get_quadrature()).WillByDefault(Return(this->quad));
      ON_CALL(*(this->sweeper.get()), status()).WillByDefault(ReturnRef(this->status));
      ON_CALL(*(this->sweeper.get()), get_status()).WillByDefault(Return(this->status));
      ON_CALL(*(this->sweeper.get()), get_initial_state()).WillByDefault(Return(this->encap));
    }
};

TEST_F(Setup, adding_coarser_level)
{
  ASSERT_THAT(controller->get_num_levels(), Eq(0));

  controller->add_sweeper(sweeper);
  EXPECT_THAT(controller->get_num_levels(), Eq(1));
}

TEST_F(Setup, a_level_must_be_added)
{
  controller->status()->t_end() = 0.2;
  controller->status()->dt() = 0.1;
  controller->status()->max_iterations() = 1;

  EXPECT_THROW(controller->setup(), std::logic_error);

  controller->add_sweeper(sweeper);

  EXPECT_CALL(*(sweeper.get()), status()).Times(1);
  EXPECT_CALL(*(sweeper.get()), setup()).Times(1);
  controller->setup();
}

TEST_F(Setup, setup_required_for_running)
{
  controller->status()->t_end() = 0.1;
  controller->status()->dt() = 0.1;
  controller->status()->max_iterations() = 1;

  controller->add_sweeper(sweeper);

  ASSERT_FALSE(controller->is_ready());
  EXPECT_THROW(controller->run(), std::logic_error);

  EXPECT_CALL(*(sweeper.get()), setup()).Times(1);
  controller->setup();

  EXPECT_TRUE(controller->is_ready());
  controller->run();
}


class Logic
  : public ::testing::Test
{
  protected:
    shared_ptr<SDC<transfer_t>>        controller;
    vector<double>                     nodes{0.0, 0.5, 1.0};
    shared_ptr<pfasst::Status<double>> status;
    shared_ptr<sweeper_t>              sweeper;
    shared_ptr<quadrature_t>           quad;
    shared_ptr<encap_t>                encap;

    virtual void SetUp()
    {
      this->encap = std::make_shared<encap_t>(vector<double>{0.0, 0.0});
      this->controller = std::make_shared<SDC<transfer_t>>();
      this->status = std::make_shared<pfasst::Status<double>>();
      this->controller->status() = status;
      this->quad = std::make_shared<quadrature_t>();
      ON_CALL(*(this->quad.get()), right_is_node()).WillByDefault(Return(true));
      ON_CALL(*(this->quad.get()), get_nodes()).WillByDefault(ReturnRef(this->nodes));
      ON_CALL(*(this->quad.get()), get_num_nodes()).WillByDefault(Return(3));

      this->sweeper = std::make_shared<sweeper_t>();
      ON_CALL(*(this->sweeper.get()), get_quadrature()).WillByDefault(Return(this->quad));
      ON_CALL(*(this->sweeper.get()), status()).WillByDefault(ReturnRef(this->status));
      ON_CALL(*(this->sweeper.get()), get_status()).WillByDefault(Return(this->status));
      ON_CALL(*(this->sweeper.get()), get_initial_state()).WillByDefault(Return(this->encap));

      this->controller->add_sweeper(this->sweeper);
    }
};

TEST_F(Logic, advance_in_time_with_sufficient_t_end)
{
  controller->status()->dt() = 0.1;
  controller->status()->time() = 1.0;
  controller->status()->step() = 1;
  controller->status()->t_end() = 2.0;

  EXPECT_CALL(*(sweeper.get()), advance(_)).Times(1);

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

  EXPECT_CALL(*(sweeper.get()), advance(_)).Times(0);

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

  EXPECT_CALL(*(sweeper.get()), advance(_)).Times(1);

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

  EXPECT_CALL(*(sweeper.get()), converged(_)).Times(1).WillRepeatedly(Return(false));

  EXPECT_FALSE(controller->advance_iteration());

  EXPECT_THAT(controller->get_status()->get_iteration(), Eq(1));
}

TEST_F(Logic, advance_iteration)
{
  controller->status()->max_iterations() = 5;
  ASSERT_THAT(controller->get_status()->get_iteration(), Eq(0));
  ASSERT_THAT(controller->get_status()->get_max_iterations(), Eq(5));

  EXPECT_CALL(*(sweeper.get()), converged(_)).Times(1);
  EXPECT_CALL(*(sweeper.get()), save()).Times(1);

  EXPECT_TRUE(controller->advance_iteration());

  EXPECT_THAT(controller->get_status()->get_iteration(), Eq(1));
}

TEST_F(Logic, single_time_step_sdc)
{
  controller->status()->max_iterations() = 3;
  controller->status()->dt() = 0.1;
  controller->status()->time() = 0.0;
  controller->status()->t_end() = 0.1;

  EXPECT_CALL(*(sweeper.get()), status()).Times(1);
  EXPECT_CALL(*(sweeper.get()), setup()).Times(1);
  controller->setup();
  Mock::VerifyAndClearExpectations(&(*(sweeper.get())));

  EXPECT_CALL(*(sweeper.get()), converged(_)).Times(4);
  EXPECT_CALL(*(sweeper.get()), save()).Times(3);

  EXPECT_CALL(*(sweeper.get()), pre_predict()).Times(1);
  EXPECT_CALL(*(sweeper.get()), predict()).Times(1);
  EXPECT_CALL(*(sweeper.get()), post_predict()).Times(1);

  EXPECT_CALL(*(sweeper.get()), pre_sweep()).Times(3);
  EXPECT_CALL(*(sweeper.get()), sweep()).Times(3);
  EXPECT_CALL(*(sweeper.get()), post_sweep()).Times(3);

  EXPECT_CALL(*(sweeper.get()), advance(_)).Times(0);
  EXPECT_CALL(*(sweeper.get()), post_step()).Times(1);

  controller->run();
  EXPECT_THAT(controller->status()->get_step(), Eq(0));
  EXPECT_THAT(controller->status()->get_iteration(), Eq(3));
}


TEST_MAIN()
