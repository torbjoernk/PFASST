#include "pfasst/comm/communicator.hpp"

#include <exception>
#include <stdexcept>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  namespace comm
  {
    template<class CommTagT>
    size_t
    Communicator<CommTagT>::get_size() const
    {
      return 0;
    }

    template<class CommTagT>
    size_t
    Communicator<CommTagT>::get_rank() const
    {
      return 0;
    }

    template<class CommTagT>
    size_t
    Communicator<CommTagT>::get_root() const
    {
      return 0;
    }

    template<class CommTagT>
    bool
    Communicator<CommTagT>::is_first() const
    {
      return false;
    }

    template<class CommTagT>
    bool
    Communicator<CommTagT>::is_last() const
    {
      return false;
    }

    template<class CommTagT>
    void
    Communicator<CommTagT>::cleanup(const bool discard)
    {
      UNUSED(discard);
    }

    template<class CommTagT>
    void
    Communicator<CommTagT>::abort(const int& err_code)
    {
      UNUSED(err_code);
      std::abort();
    }


    template<class CommTagT>
    bool
    Communicator<CommTagT>::probe(const int src_rank, const int tag)
    {
      UNUSED(src_rank); UNUSED(tag);
      throw runtime_error("not implemented: probing for incoming message");
    }


    template<class CommTagT>
    void
    Communicator<CommTagT>::send(const double* const data, const int count, const int dest_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(dest_rank); UNUSED(tag);
      throw runtime_error("not implemented: send of double");
    }

    template<class CommTagT>
    void
    Communicator<CommTagT>::send_status(const StatusDetail<double>* const data, const int count, const int dest_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(dest_rank); UNUSED(tag);
      throw runtime_error("not implemented: send of status details");
    }


    template<class CommTagT>
    void
    Communicator<CommTagT>::isend(const double* const data, const int count, const int dest_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(dest_rank); UNUSED(tag);
      throw runtime_error("not implemented: isend of double");
    }

    template<class CommTagT>
    void
    Communicator<CommTagT>::isend_status(const StatusDetail<double>* const data, const int count, const int dest_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(dest_rank); UNUSED(tag);
      throw runtime_error("not implemented: isend of status details");
    }


    template<class CommTagT>
    void
    Communicator<CommTagT>::recv(double* data, const int count, const int src_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(src_rank); UNUSED(tag);
      throw runtime_error("not implemented: recv of double");
    }

    template<class CommTagT>
    void
    Communicator<CommTagT>::recv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(src_rank); UNUSED(tag);
      throw runtime_error("not implemented: recv of status details");
    }


    template<class CommTagT>
    void
    Communicator<CommTagT>::irecv(double* data, const int count, const int src_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(src_rank); UNUSED(tag);
      throw runtime_error("not implemented: irecv of double");
    }

    template<class CommTagT>
    void
    Communicator<CommTagT>::irecv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(src_rank); UNUSED(tag);
      throw runtime_error("not implemented: irecv of status details");
    }


    template<class CommTagT>
    void
    Communicator<CommTagT>::bcast(double* data, const int count, const int root_rank)
    {
      UNUSED(data); UNUSED(count); UNUSED(root_rank);
      throw runtime_error("not implemented: bcast of double");
    }
  }  // ::pfasst::comm
}  // ::pfasst
