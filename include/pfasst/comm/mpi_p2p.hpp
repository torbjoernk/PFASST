#ifndef _PFASST__COMM__MPI_P2P_HPP_
#define _PFASST__COMM__MPI_P2P_HPP_

#ifndef WITH_MPI
  #error "You need MPI enabled for using the MPI P2P communicator"
#endif

#include <list>
#include <map>
#include <memory>
#include <type_traits>
#include <utility>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <mpi.h>
#include <leathers/pop>

#include "pfasst/comm/communicator.hpp"
#include "pfasst/controller/status.hpp"


#ifndef NON_CONST_MPI
  template<typename T>
  T* mpi_const_cast(const T* input) {
      return const_cast<T*>(input);
  }
  template<typename T>
  T& mpi_const_cast(const T& input) {
      return const_cast<T&>(input);
  }
#else
  template<typename T>
  const T* mpi_const_cast(const T* input) {
      return input;
  }
  template<typename T>
  const T& mpi_const_cast(const T& input) {
      return input;
  }
#endif


namespace pfasst
{
  namespace comm
  {
    static string error_from_code(const int err_code);
    static MPI_Status MPI_Status_factory();
    static void check_mpi_error(const int err_code);


    template<
      class CommTagT = temporal_communicator_tag
    >
    class MpiP2P
      : public Communicator<CommTagT>
    {
      protected:
        int      _size = -1;
        int      _rank = -1;
        string   _name = "";

        MPI_Comm _comm;

        vector<shared_ptr<MPI_Request>> _requests;

      public:
        using type_tag = CommTagT;

      public:
        explicit MpiP2P(MPI_Comm comm);
        MpiP2P();
        MpiP2P(const MpiP2P<CommTagT>& other) = default;
        MpiP2P(MpiP2P<CommTagT>&& other) = default;
        virtual ~MpiP2P();
        MpiP2P<CommTagT>& operator=(const MpiP2P<CommTagT>& other) = default;
        MpiP2P<CommTagT>& operator=(MpiP2P<CommTagT>&& other) = default;

        virtual size_t get_size() const override;
        virtual size_t get_rank() const override;

        virtual string get_name() const;

        virtual bool is_first() const override;
        virtual bool is_last() const override;

        virtual void cleanup(const bool discard = false) override;
        virtual void abort(const int& err_code) override;

        virtual bool probe(const int src_rank, const int tag) override;

        virtual void send(const double* const data, const int count, const int dest_rank, const int tag) override;
        virtual void send_status(const StatusDetail<double>* const data, const int count, const int dest_rank, const int tag) override;

        virtual void isend(const double* const data, const int count, const int dest_rank, const int tag) override;
        virtual void isend_status(const StatusDetail<double>* const data, const int count, const int dest_rank, const int tag) override;

        virtual void recv(double* data, const int count, const int src_rank, const int tag) override;
        virtual void recv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag) override;

        virtual void irecv(double* data, const int count, const int src_rank, const int tag) override;
        virtual void irecv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag) override;

        virtual void bcast(double* data, const int count, const int root_rank) override;
    };


    pair<shared_ptr<MpiP2P<temporal_communicator_tag>>, shared_ptr<MpiP2P<spatial_communicator_tag>>>
    split_comm(const size_t np_space);
  }  // ::pfasst::comm
}  // ::pfasst

#include "pfasst/comm/mpi_p2p_impl.hpp"

#endif  // _PFASST__COMM__MPI_P2P_HPP_
