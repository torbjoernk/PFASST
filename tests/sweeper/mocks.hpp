#include "fixtures/test_helpers.hpp"

#include <memory>
#include <vector>
using std::shared_ptr;
using std::vector;

#include "pfasst/sweeper/sweeper.hpp"


template<
  class SweeperTrait,
  typename Enabled = void
>
class SweeperMock
  : public pfasst::Sweeper<SweeperTrait, Enabled>
{
  public:
    using traits = SweeperTrait;

  public:
    MOCK_METHOD0_T(quadrature, shared_ptr<IQuadrature<typename SweeperTrait::time_t>>&());
    MOCK_CONST_METHOD0_T(get_quadrature, const shared_ptr<IQuadrature<typename SweeperTrait::time_t>>());

    MOCK_METHOD0_T(status, shared_ptr<pfasst::Status<typename SweeperTrait::time_t>>&());
    MOCK_CONST_METHOD0_T(get_status, const shared_ptr<pfasst::Status<typename SweeperTrait::time_t>>());

//     MOCK_METHOD0_T(encap_factory, shared_ptr<typename SweeperTrait::encap_t::factory_t>&());
//     MOCK_CONST_METHOD0_T(get_encap_factory, const shared_ptr<typename SweeperTrait::encap_t::factory_t>());

    MOCK_METHOD0_T(initial_state, shared_ptr<typename SweeperTrait::encap_t>&());
    MOCK_CONST_METHOD0_T(get_initial_state, const shared_ptr<typename SweeperTrait::encap_t>());
    MOCK_CONST_METHOD0_T(get_end_state, const shared_ptr<typename SweeperTrait::encap_t>());

    MOCK_METHOD0_T(setup, void());

    MOCK_METHOD0_T(pre_predict, void());
    MOCK_METHOD0_T(predict, void());
    MOCK_METHOD0_T(post_predict, void());

    MOCK_METHOD0_T(pre_sweep, void());
    MOCK_METHOD0_T(sweep, void());
    MOCK_METHOD0_T(post_sweep, void());

    MOCK_METHOD1_T(advance, void(const size_t& num_steps));
    MOCK_METHOD0_T(spread, void());
    MOCK_METHOD0_T(save, void());

    MOCK_METHOD0_T(post_step, void());

    MOCK_METHOD1_T(reevaluate, void(const bool initial_only));
    MOCK_METHOD0_T(reevaluate, void());
    MOCK_METHOD1_T(integrate, vector<shared_ptr<typename SweeperTrait::encap_t>>(const typename SweeperTrait::time_t& dt));

    MOCK_METHOD1_T(converged, bool(const bool));
    MOCK_METHOD0_T(converged, bool());
};
