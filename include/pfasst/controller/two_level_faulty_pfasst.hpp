#ifndef _PFASST__CONTROLLER__TWO_LEVEL_FAULTY_PFASST_HPP_
#define _PFASST__CONTROLLER__TWO_LEVEL_FAULTY_PFASST_HPP_

#include <memory>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <better-enums/enum.h>
#include <leathers/pop>

#include "pfasst/controller/status.hpp"
#include "pfasst/controller/two_level_pfasst.hpp"
#include "pfasst/comm/mpi_p2p.hpp"


namespace pfasst
{
  template<
    class TransferT,
    class CommT = comm::MpiP2P<>
  >
  class TwoLevelFaultyPfasst
    : public TwoLevelPfasst<TransferT, CommT>
  {
    using TagLevel = pfasst::detail::TagLevel;
    using TagModifier = pfasst::detail::TagModifier;
    using TagType = pfasst::detail::TagType;

    public:
      using transfer_t = TransferT;
      using comm_t = CommT;
      using time_t = typename transfer_t::traits::fine_time_t;

      static void init_opts();
      static void init_loggers();

    protected:
      // map: rank -> step -> iter
      map<size_t, map<size_t, size_t>> _faults;

      virtual void get_check_prev_status() override;
      virtual void recv_coarse() override;
      virtual void predictor() override;
      bool be_faulty();

      virtual int compute_tag(const TagType type,
                              const TagLevel level = TagLevel::ANY,
                              const TagModifier mod = TagModifier::UNMOD) const override;

    public:
      TwoLevelFaultyPfasst();
      TwoLevelFaultyPfasst(const TwoLevelFaultyPfasst<TransferT, CommT>& other) = default;
      TwoLevelFaultyPfasst(TwoLevelFaultyPfasst<TransferT, CommT>&& other) = default;
      virtual ~TwoLevelFaultyPfasst() = default;
      TwoLevelFaultyPfasst<TransferT, CommT>& operator=(const TwoLevelFaultyPfasst<TransferT, CommT>& other) = default;
      TwoLevelFaultyPfasst<TransferT, CommT>& operator=(TwoLevelFaultyPfasst<TransferT, CommT>&& other) = default;

      virtual void set_options() override;

      virtual void setup() override;
      virtual void run() override;

      virtual bool advance_iteration() override;
  };
}  // ::pfasst

#include "pfasst/controller/two_level_faulty_pfasst_impl.hpp"

#endif  // _PFASST__CONTROLLER__TWO_LEVEL_FAULTY_PFASST_HPP_
