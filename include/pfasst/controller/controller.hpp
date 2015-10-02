#ifndef _PFASST__CONTROLLER__INTERFACE_HPP_
#define _PFASST__CONTROLLER__INTERFACE_HPP_

#include <memory>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/comm/communicator.hpp"
#include "pfasst/controller/status.hpp"


namespace pfasst
{
  template<
    class TransferT,
    class CommT = comm::Communicator
  >
  class Controller
    : public enable_shared_from_this<Controller<TransferT, CommT>>
  {
    public:
      using transfer_t = TransferT;
      using comm_t = CommT;
      using time_t = typename transfer_t::traits::fine_time_t;

    protected:
      shared_ptr<comm_t>          _comm;

      shared_ptr<transfer_t>      _transfer;
      shared_ptr<Status<time_t>>  _status;
      bool                        _ready;
      string                      _logger_id;

      virtual void compute_num_steps();
      virtual bool& ready();

    public:
      Controller();
      Controller(const Controller<TransferT, CommT>& other) = default;
      Controller(Controller<TransferT, CommT>&& other) = default;
      virtual ~Controller() = default;
      Controller<TransferT, CommT>& operator=(const Controller<TransferT, CommT>& other) = default;
      Controller<TransferT, CommT>& operator=(Controller<TransferT, CommT>&& other) = default;

      virtual       shared_ptr<CommT>& communicator();
      virtual const shared_ptr<CommT>  get_communicator() const;

      virtual       shared_ptr<Status<typename TransferT::traits::fine_time_t>>& status();
      virtual const shared_ptr<Status<typename TransferT::traits::fine_time_t>>  get_status() const;

      virtual size_t get_num_levels() const;
      virtual bool   is_ready() const;

      virtual       void  set_logger_id(const string& logger_id);
      virtual const char* get_logger_id() const;

      template<class SweeperT>
      void add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse);

      // XXX: might not be in interface
      virtual void add_transfer(shared_ptr<TransferT> transfer);

      virtual const shared_ptr<TransferT> get_transfer() const;
      virtual       shared_ptr<TransferT> get_transfer();

      virtual void set_options();

      virtual void setup();
      virtual void run();
      virtual void post_run();

      virtual bool advance_time(const size_t& num_steps = 1);
      virtual bool advance_iteration();
  };
}  // ::pfasst

#include "pfasst/controller/controller_impl.hpp"

#endif  // _PFASST__CONTROLLER__INTERFACE_HPP_
