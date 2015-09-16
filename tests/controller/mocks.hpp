#include "fixtures/test_helpers.hpp"

#include "pfasst/controller/controller.hpp"
#include "pfasst/controller/status.hpp"
#include "pfasst/comm/communicator.hpp"


template<typename precision>
class StatusMock
  : public pfasst::Status<precision>
{
  public:
    StatusMock() = default;
    StatusMock(const StatusMock<precision>& other) = default;
    StatusMock(StatusMock<precision>&& other) = default;
    virtual ~StatusMock() = default;
    StatusMock<precision>& operator=(const StatusMock<precision>& other) = default;
    StatusMock<precision>& operator=(StatusMock<precision>&& other) = default;

    MOCK_METHOD0_T(step, size_t&());
    MOCK_CONST_METHOD0_T(get_step, size_t());

    MOCK_METHOD0_T(iteration, size_t&());
    MOCK_CONST_METHOD0_T(get_iteration, size_t());

    MOCK_METHOD0_T(time, precision&());
    MOCK_CONST_METHOD0_T(get_time, precision());

    MOCK_METHOD0_T(dt, precision&());
    MOCK_CONST_METHOD0_T(get_dt, precision());

    MOCK_METHOD1_T(set_primary_state, void(const pfasst::PrimaryState& state));
    MOCK_CONST_METHOD0_T(get_primary_state, pfasst::PrimaryState());

    MOCK_METHOD1_T(set_secondary_state, void(const pfasst::SecondaryState& state));
    MOCK_CONST_METHOD0_T(get_secondary_state, pfasst::SecondaryState());

    MOCK_METHOD0_T(residual, precision&());
    MOCK_CONST_METHOD0_T(get_residual, precision());

    MOCK_CONST_METHOD1_T(log, void(el::base::type::ostream_t& os));
};
