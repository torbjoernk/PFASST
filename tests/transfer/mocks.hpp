#include "fixtures/test_helpers.hpp"

#include <memory>
using namespace std;

#include "pfasst/transfer/transfer.hpp"


template<
  class TransferTraits,
  typename Enabled = void
>
class TransferMock
  : public pfasst::Transfer<TransferTraits, Enabled>
{
  public:
    using traits = TransferTraits;

  public:
    MOCK_METHOD2_T(interpolate_initial, void(const shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse,
                                             shared_ptr<typename TransferTraits::fine_sweeper_t> fine));
    MOCK_METHOD3_T(interpolate, void(const shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse,
                                     shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                                     const bool initial));
    MOCK_METHOD2_T(interpolate_data, void(const shared_ptr<typename TransferTraits::coarse_encap_t> coarse,
                                          shared_ptr<typename TransferTraits::fine_encap_t> fine));

    MOCK_METHOD2_T(restrict_initial, void(const shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                                    shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse));
    MOCK_METHOD3_T(restrict, void(const shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                                  shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse,
                                  const bool initial));
    MOCK_METHOD2_T(restrict_data, void(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                       shared_ptr<typename TransferTraits::coarse_encap_t> coarse));

    MOCK_METHOD3_T(fas, void(const typename TransferTraits::fine_time_t& dt,
                             const shared_ptr<typename TransferTraits::fine_sweeper_t> fine,
                             shared_ptr<typename TransferTraits::coarse_sweeper_t> coarse));
};
