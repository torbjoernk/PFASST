#include "pfasst/comm/mpi_p2p.hpp"

#include <exception>
#include <string>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"


MAKE_LOGGABLE(MPI_Status, mpi_status, os)
{
  if (   mpi_status.MPI_TAG == MPI_ANY_TAG
         && mpi_status.MPI_SOURCE == MPI_ANY_SOURCE
         && mpi_status.MPI_ERROR == MPI_SUCCESS) {
    os << "MPI_Status(empty)";
  } else {
    char err_str[MPI_MAX_ERROR_STRING];
    int err_len = 0;
    int err = MPI_Error_string(mpi_status.MPI_ERROR, err_str, &err_len);
    pfasst::comm::check_mpi_error(err);
    os << "MPI_Status(source=" << to_string(mpi_status.MPI_SOURCE) << ", "
       << "tag=" << to_string(mpi_status.MPI_TAG) << ", "
       << "error=" << string(err_str, err_len) << ")";
  }
  return os;
}


namespace pfasst
{
  namespace comm
  {
    string error_from_code(const int err_code)
    {
      char err_str[MPI_MAX_ERROR_STRING];
      int err_len = 0;
      int err = MPI_Error_string(err_code, err_str, &err_len);
      check_mpi_error(err);
      return string(err_str, err_len) + " (code=" + to_string(err_code) + ")";
    }


    MPI_Status MPI_Status_factory()
    {
      MPI_Status stat;
      stat.MPI_ERROR = MPI_SUCCESS;
      stat.MPI_SOURCE = MPI_ANY_SOURCE;
      stat.MPI_TAG = MPI_ANY_TAG;
      return stat;
    }

    void check_mpi_error(const int err_code)
    {
      if (err_code != MPI_SUCCESS) {
        string err_msg = error_from_code(err_code);
        ML_CLOG(ERROR, "COMM_P2P", "MPI encountered an error: " << err_msg);
        throw runtime_error("MPI encountered an error: " + err_msg);
      }
    }


    template<class CommTagT>
    MpiP2P<CommTagT>::MpiP2P()
      : MpiP2P(MPI_COMM_WORLD)
    {}

    template<class CommTagT>
    MpiP2P<CommTagT>::MpiP2P(MPI_Comm comm)
      :   _comm(comm)
        , _requests(0)
    {
      log::add_custom_logger("COMM_P2P");

      // get communicator's size and processors rank
      MPI_Comm_size(this->_comm, &(this->_size));
      MPI_Comm_rank(this->_comm, &(this->_rank));

      // get communicator's name (if available)
      int len = -1;
      char buff[MPI_MAX_OBJECT_NAME];
      int err = MPI_Comm_get_name(this->_comm, buff, &len);
      check_mpi_error(err);
      if (len > 0) {
        this->_name = string(buff, len);
      }
    }

    template<class CommTagT>
    MpiP2P<CommTagT>::~MpiP2P()
    {
      this->cleanup(true);
    }

    template<class CommTagT>
    size_t
    MpiP2P<CommTagT>::get_size() const
    {
      assert(this->_size > 0);
      return this->_size;
    }

    template<class CommTagT>
    size_t
    MpiP2P<CommTagT>::get_rank() const
    {
      assert(this->_rank >= 0);
      return this->_rank;
    }

    template<class CommTagT>
    string
    MpiP2P<CommTagT>::get_name() const
    {
      return this->_name;
    }

    template<class CommTagT>
    bool
    MpiP2P<CommTagT>::is_first() const
    {
      return (this->get_rank() == this->get_root());
    }

    template<class CommTagT>
    bool
    MpiP2P<CommTagT>::is_last() const
    {
      return (this->get_rank() == this->get_size() - 1);
    }

    template<class CommTagT>
    void
    MpiP2P<CommTagT>::cleanup(const bool discard)
    {
      ML_CLOG(DEBUG, "COMM_P2P",
              "cleaning up " << this->_requests.size() << " dangling request handlers");
      int err = -1;

      for (auto& req : this->_requests) {
        MPI_Status stat = MPI_Status_factory();
        err = MPI_Wait(req.get(), &stat);
        check_mpi_error(err);
        if (*(req.get()) == MPI_REQUEST_NULL) {
          req.reset();
        }
      }

      this->_requests.erase(remove(this->_requests.begin(), this->_requests.end(), nullptr),
                            this->_requests.end());

      ML_CLOG(DEBUG, "COMM_P2P", "done");
    }

    template<class CommTagT>
    void
    MpiP2P<CommTagT>::abort(const int& err_code)
    {
      MPI_Abort(this->_comm, err_code);
    }


    template<class CommTagT>
    bool
    MpiP2P<CommTagT>::probe(const int src_rank, const int tag)
    {
      ML_CLOG(DEBUG, "COMM_P2P",
              "probing for incomming message from " << src_rank << " with tag=" << tag);
      MPI_Status stat = MPI_Status_factory();
      int flag = (int)false;
      int err = MPI_Iprobe(src_rank, tag, this->_comm, &flag, &stat);
      check_mpi_error(err);
      ML_CLOG(DEBUG, "COMM_P2P", "probed: " << stat);
      return (bool)flag;
    }


    template<class CommTagT>
    void
    MpiP2P<CommTagT>::send(const double* const data, const int count, const int dest_rank, const int tag)
    {
      ML_CLOG(DEBUG, "COMM_P2P",
              "sending " << count << " " << ((count == 1) ? "double" : "doubles")
              << " with tag=" << tag << " to " << dest_rank);

      int err = MPI_Send(mpi_const_cast<void>(data), count, MPI_DOUBLE, dest_rank, tag, this->_comm);
      check_mpi_error(err);
    }

    template<class CommTagT>
    void
    MpiP2P<CommTagT>::send_status(const StatusDetail<double>* const data, const int count,
                             const int dest_rank, const int tag)
    {
      assert(pfasst::status_data_type != MPI_DATATYPE_NULL);

      ML_CLOG(DEBUG, "COMM_P2P",
              "sending " << count << " " << ((count == 1) ? "Status" : "Stati")
              << " with tag=" << tag << " to " << dest_rank);

      int err = MPI_Send(mpi_const_cast<void>(data), count, status_data_type, dest_rank, tag,
                         this->_comm);
      check_mpi_error(err);
    }


    template<class CommTagT>
    void
    MpiP2P<CommTagT>::isend(const double* const data, const int count, const int dest_rank, const int tag)
    {
      ML_CLOG(DEBUG, "COMM_P2P",
              "non-blocking send of " << count << " " << ((count == 1) ? "double" : "doubles")
              << " with tag=" << tag << " to " << dest_rank);

      this->_requests.emplace_back(new MPI_Request(MPI_REQUEST_NULL));

      int err = MPI_Isend(mpi_const_cast<void>(data), count, MPI_DOUBLE, dest_rank, tag,
                          this->_comm, this->_requests.back().get());
      check_mpi_error(err);
    }

    template<class CommTagT>
    void
    MpiP2P<CommTagT>::isend_status(const StatusDetail<double>* const data, const int count,
                              const int dest_rank, const int tag)
    {
      assert(pfasst::status_data_type != MPI_DATATYPE_NULL);

      ML_CLOG(DEBUG, "COMM_P2P",
              "non-blocking send of " << count << " " << ((count == 1) ? "Status" : "Stati")
              << " with tag=" << tag << " to " << dest_rank);

      this->_requests.emplace_back(new MPI_Request(MPI_REQUEST_NULL));

      int err = MPI_Isend(mpi_const_cast<void>(data), count, status_data_type, dest_rank, tag,
                          this->_comm, this->_requests.back().get());
      check_mpi_error(err);
    }


    template<class CommTagT>
    void
    MpiP2P<CommTagT>::recv(double* data, const int count, const int dest_rank, const int tag)
    {
      auto stat = MPI_Status_factory();
      ML_CLOG(DEBUG, "COMM_P2P",
              "receiving " << count << " " << ((count == 1) ? "double" : "doubles")
              << " with tag=" << tag << " from " << dest_rank);

      int err = MPI_Recv(data, count, MPI_DOUBLE, dest_rank, tag, this->_comm, &stat);
      check_mpi_error(err);
    }

    template<class CommTagT>
    void
    MpiP2P<CommTagT>::recv_status(StatusDetail<double>* data, const int count, const int dest_rank, const int tag)
    {
      assert(pfasst::status_data_type != MPI_DATATYPE_NULL);

      auto stat = MPI_Status_factory();
      ML_CLOG(DEBUG, "COMM_P2P",
              "receiving " << count << " " << ((count == 1) ? "Status" : "Stati")
              << " with tag=" << tag << " from " << dest_rank);

      int err = MPI_Recv(data, count, pfasst::status_data_type, dest_rank, tag, this->_comm, &stat);
      check_mpi_error(err);
    }


    template<class CommTagT>
    void
    MpiP2P<CommTagT>::irecv(double* data, const int count, const int src_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(src_rank); UNUSED(tag);
      throw runtime_error("don't use non-blocking receive on data");
    }

    template<class CommTagT>
    void
    MpiP2P<CommTagT>::irecv_status(StatusDetail<double>* data, const int count, const int src_rank, const int tag)
    {
      UNUSED(data); UNUSED(count); UNUSED(src_rank); UNUSED(tag);
      throw runtime_error("don't use non-blocking receive on status data");
    }


    template<class CommTagT>
    void
    MpiP2P<CommTagT>::bcast(double* data, const int count, const int root_rank)
    {
      ML_CLOG(DEBUG, "COMM_P2P",
              "broadcasting " << count << " " << ((count == 1) ? "double" : "doubles")
              << " from root " << root_rank);

      int err = MPI_Bcast(data, count, MPI_DOUBLE, root_rank, mpi_const_cast<MPI_Comm>(this->_comm));
      check_mpi_error(err);
    }


    pair<shared_ptr<MpiP2P<temporal_communicator_tag>>, shared_ptr<MpiP2P<spatial_communicator_tag>>>
    split_comm(const size_t np_space)
    {
      auto comm_world = make_shared<MpiP2P<temporal_communicator_tag>>(MPI_COMM_WORLD);
      int err = (int)false;

      // split into time comms
      MPI_Comm COMM_TIME;
      const int color_time = comm_world->get_rank() / (comm_world->get_size() / np_space);
      err = MPI_Comm_split(MPI_COMM_WORLD, color_time, comm_world->get_rank(), &COMM_TIME);
      pfasst::comm::check_mpi_error(err);
      auto comm_time = make_shared<MpiP2P<temporal_communicator_tag>>(COMM_TIME);

      // split into spatial comms
      MPI_Comm COMM_SPACE;
      const int color_space = comm_world->get_rank() % comm_time->get_size();
      err = MPI_Comm_split(MPI_COMM_WORLD, color_space, comm_world->get_rank(), &COMM_SPACE);
      pfasst::comm::check_mpi_error(err);
      auto comm_space = make_shared<MpiP2P<spatial_communicator_tag>>(COMM_SPACE);

      return {comm_time, comm_space};
    }
  } // ::pfasst::comm
}  // ::pfasst
