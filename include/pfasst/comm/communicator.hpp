#ifndef _PFASST__COMM__COMMUNICATOR_HPP_
#define _PFASST__COMM__COMMUNICATOR_HPP_

#include <memory>
#include <type_traits>
using namespace std;

#include "pfasst/controller/status.hpp"


namespace pfasst
{
  namespace comm
  {
    struct communicator_tag
    {};

    struct temporal_communicator_tag
      : public communicator_tag
    {};

    struct spatial_communicator_tag
      : public communicator_tag
    {};


    template<
      class CommTagT = temporal_communicator_tag
    >
    class Communicator
      : public enable_shared_from_this<Communicator<CommTagT>>
    {
      static_assert(is_base_of<communicator_tag, CommTagT>::value,
                    "Communicator must either be for temporal or spatial domain.");

      public:
        using type_tag = CommTagT;

      public:
        Communicator() = default;
        Communicator(const Communicator<CommTagT>& other) = default;
        Communicator(Communicator<CommTagT>&& other) = default;
        virtual ~Communicator() = default;
        Communicator<CommTagT>& operator=(const Communicator<CommTagT>& other) = default;
        Communicator<CommTagT>& operator=(Communicator<CommTagT>&& other) = default;

        virtual size_t get_size() const;
        virtual size_t get_rank() const;
        virtual size_t get_root() const;

        virtual bool is_first() const;
        virtual bool is_last() const;

        virtual void cleanup(const bool discard = false);
        virtual void abort(const int& err_code);

        virtual bool probe(const int src_rank, const int tag);

        virtual void send(const double* const data, const int count, const int dest_rank, const int tag);
        virtual void send_status(const StatusDetail<double>* const data, const int count, const int dest_rank, const int tag);

        virtual void isend(const double* const data, const int count, const int dest_rank, const int tag);
        virtual void isend_status(const StatusDetail<double>* const data, const int count, const int dest_rank, const int tag);

        virtual void recv(double* data, const int count, const int src_rank, const int tag);
        virtual void recv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag);

        virtual void irecv(double* data, const int count, const int src_rank, const int tag);
        virtual void irecv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag);

        virtual void bcast(double* data, const int count, const int root_rank);
    };
  }  // ::pfasst::comm
}  // ::pfasst

#include "pfasst/comm/communicator_impl.hpp"

#endif  // _PFASST__COMM__COMMUNICATOR_HPP_
