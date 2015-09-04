#include "fixtures/test_helpers.hpp"

#include <leathers/push>
#include <leathers/all>
#include <mpi.h>
#include <leathers/pop>

#include <pfasst/comm/mpi_p2p.hpp>
using pfasst::comm::MpiP2P;

INSTANTIATE_TYPED_TEST_CASE_P(MPI2P2, Concepts, ::testing::Types<MpiP2P>);


class StateAwareness
  : public ::testing::Test
{
  protected:
    MpiP2P comm;
};

TEST_F(StateAwareness, communicator_size)
{
  int size = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  EXPECT_THAT(comm.get_size(), Eq(size));
}

TEST_F(StateAwareness, rank_in_communicator)
{
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  EXPECT_THAT(comm.get_rank(), Eq(rank));
}

TEST_F(StateAwareness, root_rank)
{
  EXPECT_THAT(comm.get_root(), Eq(0));
}

TEST_F(StateAwareness, being_first_or_last)
{
  if (comm.get_rank() == comm.get_root()) {
    EXPECT_TRUE(comm.is_first());
    EXPECT_FALSE(comm.is_last());
  } else if (comm.get_rank() == comm.get_size() - 1) {
    EXPECT_TRUE(comm.is_last());
    EXPECT_FALSE(comm.is_first());
  } else {
    EXPECT_FALSE(comm.is_first());
    EXPECT_FALSE(comm.is_last());
  }
}



class MessagePassing
  : public ::testing::Test
{
  protected:
    MpiP2P comm;
};

TEST_F(MessagePassing, blocking_send_recv_double)
{
  double value = 0.0;
  if (comm.get_rank() % 2 == 0) {
    value = 42.21;
    comm.send(&value, 1, comm.get_rank() + 1, 0);
  } else if (comm.get_rank() % 2 == 1) {
    comm.recv(&value, 1, comm.get_rank() - 1, 0);
  }
  EXPECT_THAT(value, DoubleEq(42.21));
}

TEST_F(MessagePassing, non_blocking_send_recv_double)
{
  double value = 0.0;
  if (comm.get_rank() % 2 == 0) {
    value = 42.21;
    comm.isend(&value, 1, comm.get_rank() + 1, 0);
  } else if (comm.get_rank() % 2 == 1) {
    comm.irecv(&value, 1, comm.get_rank() - 1, 0);
  }

  // this is required to wait for pending requests
  //  otherwise the assert will fail as the pending requests will get cleaned at MPI_Finalize way
  //  after the check
  comm.cleanup();

  EXPECT_THAT(value, DoubleEq(42.21));
}

TEST_F(MessagePassing, broadcast_double)
{
  double value = 0.0;
  if (comm.get_rank() == 0) {
    value = 42.21;
  }
  comm.bcast(&value, 1, 0);

  EXPECT_THAT(value, DoubleEq(42.21));
}


TEST_MAIN()
