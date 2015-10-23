#include "fixtures/test_helpers.hpp"
using ::testing::Eq;
using ::testing::Pointwise;

#include <pfasst/contrib/spectral_transfer.hpp>
using pfasst::contrib::SpectralTransfer;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
using encap_traits_t = pfasst::encap::vector_encap_traits<double, double, 3>;
using encap_t = pfasst::encap::Encapsulation<encap_traits_t>;

#include "pfasst/sweeper/sweeper.hpp"
using sweeper_t = pfasst::Sweeper<pfasst::sweeper_traits<encap_traits_t>>;

using transfer_t = SpectralTransfer<pfasst::transfer_traits<sweeper_t, sweeper_t, 2>>;

using Spectral3DTransferTypes = ::testing::Types<SpectralTransfer<pfasst::transfer_traits<sweeper_t, sweeper_t, 2>>>;
INSTANTIATE_TYPED_TEST_CASE_P(Spectral3DTransfer, Concepts, Spectral3DTransferTypes);


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
      this->coarse_encap = make_shared<encap_t>(vector<double>(8, 1.0));
      this->fine_encap = make_shared<encap_t>(64);

      ASSERT_THAT(this->coarse_encap->get_dimwise_num_dofs(), Pointwise(Eq(), array<int, 3>{{2, 2, 2}}));
      ASSERT_THAT(this->fine_encap->get_dimwise_num_dofs(), Pointwise(Eq(), array<int, 3>{{4, 4, 4}}));
    }
};

TEST_F(Interpolation, interpolate_constant)
{
  LOG(INFO) << "coarse ("<<coarse_encap->get_dimwise_num_dofs()<<"): " << coarse_encap->get_data();
  LOG(INFO) << "fine ("<<fine_encap->get_dimwise_num_dofs()<<"): " << fine_encap->get_data();
  transfer.interpolate_data(this->coarse_encap, this->fine_encap);
  EXPECT_THAT(this->fine_encap->data(), Pointwise(Eq(), vector<double>(64, 1.0)));
}


TEST_MAIN()
