#include "fixtures/test_helpers.hpp"

#include <vector>
using namespace std;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/encapsulation.hpp>
#include <pfasst/encap/vector.hpp>
using encap_traits_t = pfasst::encap::vector_encap_traits<double, double, 1>;
using encap_t = pfasst::encap::Encapsulation<encap_traits_t>;
using encap_factory_t = pfasst::encap::EncapsulationFactory<encap_traits_t>;

using FactoryTypes = ::testing::Types<encap_factory_t>;
INSTANTIATE_TYPED_TEST_CASE_P(VectorEncapFactory, Concepts, FactoryTypes);


TEST(VectorEncapFactory, set_size_after_initialization)
{
  encap_factory_t factory;
  EXPECT_THAT(factory.size(), Eq(0));

  factory.set_size(3);
  EXPECT_THAT(factory.size(), Eq(3));
}

TEST(VectorEncapFactory, takes_fixed_size)
{
  encap_factory_t factory(3);
  EXPECT_THAT(factory.size(), Eq(3));
}

TEST(VectorEncapFactory, produces_encapsulated_vectors)
{
  encap_factory_t factory(3);
  auto encap = factory.create();
  EXPECT_THAT(encap, NotNull());
  EXPECT_THAT(encap->get_data(), SizeIs(3));
}

TEST_MAIN()
