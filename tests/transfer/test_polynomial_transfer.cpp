#include "fixtures/test_helpers.hpp"

#include <pfasst/transfer/polynomial.hpp>
using pfasst::PolynomialTransfer;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
using encap_traits_t = pfasst::encap::vector_encap_traits<double, double, 1>;
using encap_t = pfasst::encap::Encapsulation<encap_traits_t>;

#include "pfasst/sweeper/sweeper.hpp"
using sweeper_t = pfasst::Sweeper<pfasst::sweeper_traits<encap_traits_t>>;

using transfer_traits_t = pfasst::transfer_traits<sweeper_t, sweeper_t, 2>;
using transfer_t = PolynomialTransfer<transfer_traits_t>;

using PolynomialTransferTypes = ::testing::Types<transfer_t>;
INSTANTIATE_TYPED_TEST_CASE_P(PolynomialTransfer, Concepts, PolynomialTransferTypes);



TEST_MAIN()
