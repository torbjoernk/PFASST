#include "fixtures/test_helpers.hpp"

#include <pfasst/contrib/spectral_transfer.hpp>
using pfasst::contrib::SpectralTransfer;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
using encap_traits_t = pfasst::encap::vector_encap_traits<double, double, 1>;
using encap_t = pfasst::encap::Encapsulation<encap_traits_t>;

#include "pfasst/sweeper/sweeper.hpp"
using sweeper_t = pfasst::Sweeper<pfasst::sweeper_traits<encap_traits_t>>;

using transfer_t = SpectralTransfer<pfasst::transfer_traits<sweeper_t, sweeper_t, 2>>;

using Spectral1DTransferTypes = ::testing::Types<SpectralTransfer<pfasst::transfer_traits<sweeper_t, sweeper_t, 2>>>;
INSTANTIATE_TYPED_TEST_CASE_P(Spectral1DTransfer, Concepts, Spectral1DTransferTypes);


class Interpolation
  : public ::testing::Test
{
  protected:
    transfer_t            transfer;
    shared_ptr<sweeper_t> coarse_sweeper;
    shared_ptr<sweeper_t> fine_sweeper;
    shared_ptr<encap_t>   coarse_encap;
    shared_ptr<encap_t>   fine_encap;

    virtual void SetUp()
    {
      this->coarse_encap = make_shared<encap_t>(vector<double>(3, 1.0));
      this->fine_encap = make_shared<encap_t>(6);

      ASSERT_THAT(this->coarse_encap->get_dimwise_num_dofs(), Pointwise(Eq(), array<int, 1>{{3}}));
      ASSERT_THAT(this->fine_encap->get_dimwise_num_dofs(), Pointwise(Eq(), array<int, 1>{{6}}));
    }
};

TEST_F(Interpolation, interpolate_constant)
{
  transfer.interpolate_data(this->coarse_encap, this->fine_encap);
  EXPECT_THAT(this->fine_encap->data(), Pointwise(Eq(), vector<double>(6, 1.0)));
}


TEST_MAIN()
