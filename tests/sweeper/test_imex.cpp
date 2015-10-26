#include "fixtures/test_helpers.hpp"
using ::testing::Each;
using ::testing::Eq;
using ::testing::IsEmpty;
using ::testing::IsNull;
using ::testing::NiceMock;
using ::testing::Not;
using ::testing::NotNull;
using ::testing::Pointwise;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::SizeIs;

#include <memory>
#include <type_traits>
#include <vector>
using std::shared_ptr;
using std::vector;

#include <pfasst/sweeper/imex.hpp>
using pfasst::IMEX;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
using encap_traits_t = pfasst::encap::vector_encap_traits<double, double, 1>;
using encap_t = pfasst::encap::Encapsulation<encap_traits_t>;
using sweeper_t = IMEX<pfasst::sweeper_traits<encap_traits_t>>;
static_assert(std::is_same<encap_t, typename sweeper_t::traits::encap_t>::value, "");

#include "quadrature/mocks.hpp"
#include "controller/mocks.hpp"

using SweeperTypes = ::testing::Types<sweeper_t>;
INSTANTIATE_TYPED_TEST_CASE_P(IMEX, Concepts, SweeperTypes);


class Setup
  : public ::testing::Test
{
  protected:
    sweeper_t                                    sweeper;
    vector<double>                               nodes{0.0, 0.5, 1.0};
    shared_ptr<NiceMock<QuadratureMock<double>>> quadrature;
    shared_ptr<NiceMock<StatusMock<double>>>     status;

    virtual void SetUp()
    {
      this->status = std::make_shared<NiceMock<StatusMock<double>>>();
      this->quadrature = std::make_shared<NiceMock<QuadratureMock<double>>>();
      ON_CALL(*(quadrature.get()), get_num_nodes()).WillByDefault(Return(nodes.size()));
      ON_CALL(*(quadrature.get()), get_nodes()).WillByDefault(ReturnRef(nodes));
    }
};

TEST_F(Setup, quadrature_is_required_for_setup)
{
  ASSERT_THAT(sweeper.quadrature(), IsNull());
  ASSERT_THAT(sweeper.get_quadrature(), IsNull());

  ASSERT_THAT(sweeper.status(), IsNull());
  ASSERT_THAT(sweeper.get_status(), IsNull());

  EXPECT_THROW(sweeper.setup(), std::runtime_error);

  sweeper.quadrature() = quadrature;
  EXPECT_THAT(sweeper.quadrature(), NotNull());
  EXPECT_THAT(sweeper.get_quadrature(), NotNull());

  sweeper.status() = status;
  EXPECT_THAT(sweeper.status(), NotNull());
  EXPECT_THAT(sweeper.get_status(), NotNull());

  sweeper.setup();
}

TEST_F(Setup, state_data_initialized_after_setup)
{
  EXPECT_THROW(sweeper.get_initial_state(), std::runtime_error);
  EXPECT_THROW(sweeper.initial_state(), std::runtime_error);
  EXPECT_THAT(sweeper.get_end_state(), IsNull());
  EXPECT_THAT(sweeper.get_states(), IsEmpty());
  EXPECT_THAT(sweeper.get_previous_states(), IsEmpty());
  EXPECT_THAT(sweeper.tau(), IsEmpty());
  EXPECT_THAT(sweeper.get_residuals(), IsEmpty());

  sweeper.status() = status;
  sweeper.quadrature() = quadrature;
  auto num_nodes = quadrature->get_num_nodes();
  sweeper.setup();

  EXPECT_THAT(sweeper.get_initial_state(), NotNull());
  EXPECT_THAT(sweeper.get_end_state(), NotNull());

  EXPECT_THAT(sweeper.get_states(), SizeIs(num_nodes + 1));
  EXPECT_THAT(sweeper.get_states(), Each(NotNull()));
  EXPECT_THAT(sweeper.get_states(), Not(MutuallyEqual()));

  EXPECT_THAT(sweeper.get_previous_states(), SizeIs(num_nodes + 1));
  EXPECT_THAT(sweeper.get_previous_states(), Each(NotNull()));
  EXPECT_THAT(sweeper.get_previous_states(), Not(MutuallyEqual()));

  EXPECT_THAT(sweeper.get_tau(), SizeIs(num_nodes + 1));
  EXPECT_THAT(sweeper.get_tau(), Each(NotNull()));
  EXPECT_THAT(sweeper.get_tau(), Not(MutuallyEqual()));

  EXPECT_THAT(sweeper.get_residuals(), SizeIs(num_nodes + 1));
  EXPECT_THAT(sweeper.get_residuals(), Each(NotNull()));
  EXPECT_THAT(sweeper.get_residuals(), Not(MutuallyEqual()));
}


class DataAccess
  : public ::testing::Test
{
  protected:
    sweeper_t                                    sweeper;
    vector<double>                               nodes{0.0, 0.5, 1.0};
    shared_ptr<NiceMock<QuadratureMock<double>>> quadrature;
    shared_ptr<encap_t>                          encap;
    shared_ptr<NiceMock<StatusMock<double>>>     status;

    virtual void SetUp()
    {
      this->status = std::make_shared<NiceMock<StatusMock<double>>>();
      this->encap = std::make_shared<encap_t>(vector<double>{1.0, 2.0, 3.0});
      this->quadrature = std::make_shared<NiceMock<QuadratureMock<double>>>();
      ON_CALL(*(quadrature.get()), get_num_nodes()).WillByDefault(Return(nodes.size()));
      ON_CALL(*(quadrature.get()), get_nodes()).WillByDefault(ReturnRef(nodes));
      sweeper.encap_factory()->set_size(3);
      sweeper.quadrature() = quadrature;
      sweeper.status() = status;
      sweeper.setup();
    }
};

TEST_F(DataAccess, initial_state_for_modification)
{
  sweeper.initial_state() = encap;
  EXPECT_THAT(sweeper.get_initial_state()->data(), Pointwise(Eq(), encap->data()));

  *(sweeper.get_initial_state().get()) = vector<double>{1.0, 1.0, 1.0};
  EXPECT_THAT(sweeper.get_initial_state()->data(), Each(1.0));
}

TEST_F(DataAccess, tau_for_modification)
{
  sweeper.tau() = vector<shared_ptr<encap_t>>{encap, encap, encap};
  EXPECT_THAT(sweeper.get_tau(), Each(Eq(encap)));
}

TEST_F(DataAccess, states_after_spreading_initial_state)
{
  sweeper.initial_state() = encap;
  sweeper.spread();
  EXPECT_THAT(sweeper.get_initial_state(), Eq(encap));
  EXPECT_THAT(sweeper.get_initial_state()->data(), Pointwise(Eq(), encap->data()));
  EXPECT_THAT(sweeper.get_states(), Not(MutuallyEqual()));
  for (auto state : sweeper.get_states()) {
    EXPECT_THAT(state->data(), Pointwise(Eq(), encap->data()));
  }
}

TEST_F(DataAccess, previous_states_after_spreading_initial_state_and_saving)
{
  sweeper.initial_state() = encap;
  sweeper.spread();
  sweeper.save();
  EXPECT_THAT(sweeper.get_initial_state(), Eq(encap));
  EXPECT_THAT(sweeper.get_initial_state()->data(), Pointwise(Eq(), encap->data()));
  EXPECT_THAT(sweeper.get_previous_states(), Not(MutuallyEqual()));
  for (auto state : sweeper.get_previous_states()) {
    EXPECT_THAT(state->data(), Pointwise(Eq(), encap->data()));
  }
}


class Logic
  : public ::testing::Test
{
  protected:
    sweeper_t                                    sweeper;
    vector<double>                               nodes{0.0, 0.5, 1.0};
    shared_ptr<NiceMock<QuadratureMock<double>>> quadrature;
    shared_ptr<NiceMock<StatusMock<double>>>     status;
    shared_ptr<encap_t>                          encap;
    Matrix<double>                               b_mat;

    virtual void SetUp()
    {
      this->b_mat = Matrix<double>::Zero(1, 4);
      this->status = std::make_shared<NiceMock<StatusMock<double>>>();
      this->encap = std::make_shared<encap_t>(vector<double>{1.0, 2.0, 3.0});
      this->quadrature = std::make_shared<NiceMock<QuadratureMock<double>>>();
      ON_CALL(*(quadrature.get()), get_num_nodes()).WillByDefault(Return(nodes.size()));
      ON_CALL(*(quadrature.get()), get_nodes()).WillByDefault(ReturnRef(nodes));

      sweeper.encap_factory()->set_size(3);
      for(size_t n = 1; n < 4; ++n) {
        b_mat(0, n) = 1.0 / 3.0;
      }
      ON_CALL(*(quadrature.get()), get_b_mat()).WillByDefault(ReturnRef(b_mat));
      sweeper.quadrature() = quadrature;
      ON_CALL(*(status.get()), get_dt()).WillByDefault(Return(1.0));
      sweeper.status() = status;
      sweeper.setup();

    }
};

TEST_F(Logic, pre_predict_copies_initial_state_if_left_is_node)
{
  ON_CALL(*(quadrature.get()), left_is_node()).WillByDefault(Return(true));

  sweeper.initial_state() = encap;
  sweeper.pre_predict();

  EXPECT_THAT(sweeper.get_states().front()->data(), Pointwise(Eq(), encap->data()));
}

TEST_F(Logic, DISABLED_pre_predict_failes_initial_state_if_left_is_node)
{
  ON_CALL(*(quadrature.get()), left_is_node()).WillByDefault(Return(false));

  sweeper.initial_state() = encap;
  sweeper.pre_predict();
}

TEST_F(Logic, post_predict_finishes_end_state_if_right_is_node)
{
  ON_CALL(*(quadrature.get()), right_is_node()).WillByDefault(Return(true));

  sweeper.initial_state() = encap;
  sweeper.spread();

  EXPECT_THAT(sweeper.get_end_state(), NotNull());
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Not(Eq()), encap->data()));

  sweeper.post_predict();
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Eq(), encap->data()));
}

TEST_F(Logic, post_predict_finishes_end_state_if_right_is_not_node)
{
  ON_CALL(*(quadrature.get()), right_is_node()).WillByDefault(Return(false));

  sweeper.initial_state() = encap;
  sweeper.spread();

  EXPECT_THAT(sweeper.get_end_state(), NotNull());
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Not(Eq()), encap->data()));

  sweeper.post_predict();
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Eq(), encap->data()));
}

TEST_F(Logic, post_sweep_finishes_end_state_if_right_is_node)
{
  ON_CALL(*(quadrature.get()), right_is_node()).WillByDefault(Return(true));

  sweeper.initial_state() = encap;
  sweeper.spread();

  EXPECT_THAT(sweeper.get_end_state(), NotNull());
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Not(Eq()), encap->data()));

  sweeper.post_sweep();
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Eq(), encap->data()));
}

TEST_F(Logic, post_sweep_finishes_end_state_if_right_is_not_node)
{
  ON_CALL(*(quadrature.get()), right_is_node()).WillByDefault(Return(false));

  sweeper.initial_state() = encap;
  sweeper.spread();

  EXPECT_THAT(sweeper.get_end_state(), NotNull());
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Not(Eq()), encap->data()));

  sweeper.post_sweep();
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Eq(), encap->data()));
}


TEST_MAIN()
