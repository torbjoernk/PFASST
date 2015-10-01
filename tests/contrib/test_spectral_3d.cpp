#include "fixtures/test_helpers.hpp"

#include <pfasst/contrib/spectral_transfer.hpp>
using pfasst::contrib::SpectralTransfer;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
using VectorEncapTrait = pfasst::vector_encap_traits<double, double, 3>;

#include "pfasst/sweeper/sweeper.hpp"
using Sweeper = pfasst::Sweeper<pfasst::sweeper_traits<VectorEncapTrait>>;

using SpectralTransferTypes = ::testing::Types<SpectralTransfer<pfasst::transfer_traits<Sweeper, Sweeper, 2>>>;
INSTANTIATE_TYPED_TEST_CASE_P(Spectral3DTransfer, Concepts, SpectralTransferTypes);


class Interpolation
  : public ::testing::Test
{
  protected:
    typedef SpectralTransfer<pfasst::transfer_traits<Sweeper, Sweeper, 2>> transfer_type;
    typedef pfasst::encap::VectorEncapsulation<double, double, 3>          encap_type;

    transfer_type          transfer;
    shared_ptr<Sweeper>    coarse_sweeper;
    shared_ptr<Sweeper>    fine_sweeper;
    shared_ptr<encap_type> coarse_encap;
    shared_ptr<encap_type> fine_encap;

    virtual void SetUp()
    {
      this->coarse_encap = make_shared<encap_type>(vector<double>(8, 1.0));
      this->fine_encap = make_shared<encap_type>(64);

      ASSERT_THAT(this->coarse_encap->get_dimwise_num_dofs(), Pointwise(Eq(), array<int, 3>{2, 2, 2}));
      ASSERT_THAT(this->fine_encap->get_dimwise_num_dofs(), Pointwise(Eq(), array<int, 3>{4, 4, 4}));
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
