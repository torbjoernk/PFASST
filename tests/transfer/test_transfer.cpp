#include "fixtures/test_helpers.hpp"

#include <stdexcept>
using namespace std;

#include <pfasst/transfer/transfer.hpp>
using pfasst::Transfer;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
using encap_traits_t = pfasst::encap::vector_encap_traits<double, double, 1>;
using encap_t = pfasst::encap::Encapsulation<encap_traits_t>;

#include "pfasst/sweeper/sweeper.hpp"
using sweeper_t = pfasst::Sweeper<pfasst::sweeper_traits<encap_traits_t>>;

using transfer_t = Transfer<pfasst::transfer_traits<sweeper_t, sweeper_t, 2>>;

using TransferTypes = ::testing::Types<transfer_t>;
INSTANTIATE_TYPED_TEST_CASE_P(Transfer, Concepts, TransferTypes);


class Interface
  : public ::testing::Test
{
  protected:
    transfer_t            transfer;
    shared_ptr<sweeper_t> coarse_sweeper;
    shared_ptr<sweeper_t> fine_sweeper;
    shared_ptr<encap_t>   coarse_encap;
    shared_ptr<encap_t>   fine_encap;
};

TEST_F(Interface, no_implementation_of_interpolation_of_initial_value)
{
  EXPECT_THROW(transfer.interpolate_initial(coarse_sweeper, fine_sweeper), runtime_error);
}

TEST_F(Interface, no_implementation_of_interpolation)
{
  EXPECT_THROW(transfer.interpolate(coarse_sweeper, fine_sweeper), runtime_error);
  EXPECT_THROW(transfer.interpolate(coarse_sweeper, fine_sweeper, true), runtime_error);
  EXPECT_THROW(transfer.interpolate(coarse_sweeper, fine_sweeper, false), runtime_error);
}

TEST_F(Interface, no_implementation_of_interpolating_encaps)
{
  EXPECT_THROW(transfer.interpolate_data(coarse_encap, fine_encap), runtime_error);
}

TEST_F(Interface, no_implementation_of_restriction_of_initial_value)
{
  EXPECT_THROW(transfer.restrict_initial(coarse_sweeper, fine_sweeper), runtime_error);
}

TEST_F(Interface, no_implementation_of_restriction)
{
  EXPECT_THROW(transfer.restrict(coarse_sweeper, fine_sweeper), runtime_error);
  EXPECT_THROW(transfer.restrict(coarse_sweeper, fine_sweeper, true), runtime_error);
  EXPECT_THROW(transfer.restrict(coarse_sweeper, fine_sweeper, false), runtime_error);
}

TEST_F(Interface, no_implementation_of_restricting_encaps)
{
  EXPECT_THROW(transfer.restrict_data(coarse_encap, fine_encap), runtime_error);
}

TEST_F(Interface, no_implementation_of_fas_correction)
{
  EXPECT_THROW(transfer.fas(1.0, fine_sweeper, coarse_sweeper), runtime_error);
}


TEST_MAIN()
