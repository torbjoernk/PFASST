#ifndef _PFASST__CONTROLLER__STATUS_HPP_
#define _PFASST__CONTROLLER__STATUS_HPP_

#include <limits>
#include <iostream>
#include <map>
#include <memory>
#include <type_traits>
#include <utility>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <better-enums/enum.h>
#include <leathers/pop>

#ifdef WITH_MPI
  #include <leathers/push>
  #include <leathers/all>
  #include <mpi.h>
  #include <leathers/pop>
#endif

#include "pfasst/logging.hpp"


namespace pfasst
{
#ifdef WITH_MPI
  static MPI_Datatype status_data_type;
#endif

  ENUM(PrimaryState, int, 
    // general state
    CONVERGED         =  0,
    FAILED            =  1,

    // iterating states
    PREDICTING        = 10,
    ITERATING         = 20,

    // intermediate states
    INTER_ITER        = 30,

    UNKNOWN_PRIMARY   = numeric_limits<int>::max()
  )

  template<class CharT, class Traits>
  inline basic_ostream<CharT, Traits>&
  operator<<(basic_ostream<CharT, Traits>& os, const pfasst::PrimaryState& state)
  {
    os << (+state)._to_string();
    return os;
  }

  ENUM(SecondaryState, int,
    // inter-level states
    CYCLE_DOWN        =  1,
    CYCLE_UP          =  2,

    // coarse level
    PRE_ITER_COARSE   =  9,
    ITER_COARSE       = 10,
    POST_ITER_COARSE  = 11,

    // fine level
    PRE_ITER_FINE     = 19,
    ITER_FINE         = 20,
    POST_ITER_FINE    = 21,

    // intermediate states
    CONV_CHECK        = 30,

    UNKNOWN_SECONDARY = numeric_limits<int>::max()
  )

  template<class CharT, class Traits>
  inline basic_ostream<CharT, Traits>&
  operator<<(basic_ostream<CharT, Traits>& os, const pfasst::SecondaryState& state)
  {
    os << (+state)._to_string();
    return os;
  }

  // the idea for this C++11 constexpr stuff is based on http://stackoverflow.com/a/26079954/588243
  using PrimarySecondaryMapItem = pair<PrimaryState, SecondaryState>;
  constexpr PrimarySecondaryMapItem VALID_STATES_COMBINATIONS[] = {
    { PrimaryState::CONVERGED,       SecondaryState::UNKNOWN_SECONDARY },

    { PrimaryState::FAILED,          SecondaryState::UNKNOWN_SECONDARY },

    { PrimaryState::PREDICTING,      SecondaryState::UNKNOWN_SECONDARY },
    { PrimaryState::PREDICTING,      SecondaryState::CYCLE_DOWN },
    { PrimaryState::PREDICTING,      SecondaryState::PRE_ITER_COARSE },
    { PrimaryState::PREDICTING,      SecondaryState::ITER_COARSE },
    { PrimaryState::PREDICTING,      SecondaryState::POST_ITER_COARSE },
    { PrimaryState::PREDICTING,      SecondaryState::CYCLE_UP },
    { PrimaryState::PREDICTING,      SecondaryState::PRE_ITER_FINE },
    { PrimaryState::PREDICTING,      SecondaryState::ITER_FINE },
    { PrimaryState::PREDICTING,      SecondaryState::POST_ITER_FINE },

    { PrimaryState::ITERATING,       SecondaryState::UNKNOWN_SECONDARY },
    { PrimaryState::ITERATING,       SecondaryState::CYCLE_DOWN },
    { PrimaryState::ITERATING,       SecondaryState::PRE_ITER_COARSE },
    { PrimaryState::ITERATING,       SecondaryState::ITER_COARSE },
    { PrimaryState::ITERATING,       SecondaryState::POST_ITER_COARSE },
    { PrimaryState::ITERATING,       SecondaryState::CYCLE_UP },
    { PrimaryState::ITERATING,       SecondaryState::PRE_ITER_FINE },
    { PrimaryState::ITERATING,       SecondaryState::ITER_FINE },
    { PrimaryState::ITERATING,       SecondaryState::POST_ITER_FINE },

    { PrimaryState::INTER_ITER,      SecondaryState::UNKNOWN_SECONDARY },
    { PrimaryState::INTER_ITER,      SecondaryState::CONV_CHECK },

    { PrimaryState::UNKNOWN_PRIMARY, SecondaryState::UNKNOWN_SECONDARY }
  };
  constexpr auto VALID_STATES_COMBINATIONS_SIZE = sizeof VALID_STATES_COMBINATIONS / sizeof VALID_STATES_COMBINATIONS[0];

  static constexpr bool validate_state_combination(PrimaryState primary, SecondaryState secondary,
                                                   int range = VALID_STATES_COMBINATIONS_SIZE) {
      return (range == 0)
               ? false
               : ((VALID_STATES_COMBINATIONS[range - 1].first == primary)
                  ? ((VALID_STATES_COMBINATIONS[range - 1].second == secondary)
                     ? true
                     : validate_state_combination(primary, secondary, range - 1)
                    )
                  : validate_state_combination(primary, secondary, range - 1)
                 );
  }
  // ... end of constexpr stuff


  template<
    typename precision
  >
  struct StatusDetail
  {
    PrimaryState   primary_state   = (+pfasst::PrimaryState::UNKNOWN_PRIMARY);
    SecondaryState secondary_state = (+pfasst::SecondaryState::UNKNOWN_SECONDARY);

    size_t    step                 = 0;
    size_t    num_steps            = 0;
    size_t    iteration            = 0;
    size_t    max_iterations       = 0;

    precision time                 = 0.0;
    precision dt                   = 0.0;
    precision t_end                = 0.0;
    precision abs_res_norm         = 0.0;
    precision rel_res_norm         = 0.0;
  };

  template<typename precision>
  inline MAKE_LOGGABLE(StatusDetail<precision>, status_detail, os)
  {
    os << LOG_FIXED << "StatusDetail("
       << "t=" << status_detail.time
       << ", dt=" << status_detail.dt
       << ", t_end=" << status_detail.t_end
       << ", step=" << status_detail.step
       << ", num_steps=" << status_detail.num_steps
       << ", iter=" << status_detail.iteration
       << ", iter_max=" << status_detail.max_iterations
       << ", 1st_state=" << (+status_detail.primary_state)._to_string()
       << ", 2nd_state=" << (+status_detail.secondary_state)._to_string()
       << ", abs_res=" << LOG_FLOAT << status_detail.abs_res_norm
       << ", rel_res=" << LOG_FLOAT << status_detail.rel_res_norm
       << ")";
    return os;
  }


  template<
    typename precision
  >
  class Status
    : public enable_shared_from_this<Status<precision>>
      , public el::Loggable
  {
    static_assert(is_arithmetic<precision>::value,
                  "precision type must be arithmetic");

    public:
      using precision_t = precision;

#ifdef WITH_MPI
      static inline void create_mpi_datatype();
      static inline void free_mpi_datatype();
#endif

      static_assert(is_standard_layout<StatusDetail<precision>>::value,
                    "Status::Detail needs to be have standard layout for MPI derived datatype");
      StatusDetail<precision> _detail;

    public:
      Status() = default;
      Status(const Status<precision>& other) = default;
      Status(Status<precision>&& other) = default;
      virtual ~Status() = default;
      Status<precision>& operator=(const Status<precision>& other) = default;
      Status<precision>& operator=(Status<precision>&& other) = default;

      virtual void clear();

      virtual size_t& step();
      virtual size_t  get_step() const;

      virtual size_t& num_steps();
      virtual size_t  get_num_steps() const;

      virtual size_t& iteration();
      virtual size_t  get_iteration() const;

      virtual size_t& max_iterations();
      virtual size_t  get_max_iterations() const;

      virtual precision& time();
      virtual precision  get_time() const;

      virtual precision& dt();
      virtual precision  get_dt() const;

      virtual precision& t_end();
      virtual precision  get_t_end() const;

      virtual void         set_primary_state(const PrimaryState& state);
      virtual PrimaryState get_primary_state() const;

      virtual void           set_secondary_state(const SecondaryState& state);
      virtual SecondaryState get_secondary_state() const;

      virtual precision& abs_res_norm();
      virtual precision  get_abs_res_norm() const;

      virtual precision& rel_res_norm();
      virtual precision  get_rel_res_norm() const;

      template<class CommT>
      bool probe(shared_ptr<CommT> comm, const int src_rank, const int tag);

      template<class CommT>
      void send(shared_ptr<CommT> comm, const int dest_rank, const int tag, const bool blocking);

      template<class CommT>
      void recv(shared_ptr<CommT> comm, const int src_rank, const int tag, const bool blocking);

      virtual vector<string> summary() const;
      virtual void log(el::base::type::ostream_t& os) const;
  };
}  // ::pfasst


#include "pfasst/controller/status_impl.hpp"

#endif  // _PFASST__CONTROLLER__STATUS_HPP_
