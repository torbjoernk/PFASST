#include "fixtures/test_helpers.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <vector>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <boost/math/constants/constants.hpp>
#include <leathers/pop>
using namespace boost::math::constants;


#include <pfasst/contrib/fft.hpp>
using pfasst::contrib::FFT;

#include <pfasst/encap/vector.hpp>
using encap_t = pfasst::encap::VectorEncapsulation<double, double, 3>;

using fft_t = FFT<encap_t>;

#include <pfasst/logging.hpp>


using FFTTypes = ::testing::Types<fft_t>;
INSTANTIATE_TYPED_TEST_CASE_P(FFT2D, Concepts, FFTTypes);

class Interface
  : public ::testing::Test
{
  protected:
    fft_t fft;
};

TEST_F(Interface, query_z_pointer_for_specific_num_dofs)
{
  complex<double>* z_ptr = fft.get_workspace({1, 1, 1})->z;
  EXPECT_THAT(z_ptr, NotNull());
}


class DiscreteFastFourierTransform
  : public ::testing::TestWithParam<shared_ptr<encap_t>>
{
  protected:
    fft_t fft;
    shared_ptr<encap_t> values;

    virtual void SetUp()
    {
      this->values = GetParam();
    }

    vector<double> two_pi_k_t(const shared_ptr<encap_t> vec, const size_t& k)
    {
      const size_t ndofs = values->get_data().size();
      vector<double> result(ndofs);

      transform(vec->get_data().cbegin(), vec->get_data().cend(),
                result.begin(),
                [k](const double& t) {
                  return cos(two_pi<double>() * k * t);
                });

      return result;
    }
};


TEST_P(DiscreteFastFourierTransform, backward_transform)
{
  const size_t ndofs = values->get_data().size();
  const size_t dim_ndofs = cbrt(ndofs);
  for (size_t k = 0; k < dim_ndofs; ++k) {
    const double precision = k * ndofs * 1.3 * numeric_limits<double>::epsilon();

    shared_ptr<encap_t> test_values = make_shared<encap_t>(two_pi_k_t(values, k));

    fft.forward(test_values);

    shared_ptr<encap_t> backward = make_shared<encap_t>();
    backward->data().resize(ndofs);
    backward->zero();
    shared_ptr<encap_t> expected = make_shared<encap_t>(backward->get_data());
    expected->scaled_add(ndofs, test_values);

    fft.backward(backward);

    for (size_t i = 0; i < ndofs; ++i) {
      EXPECT_THAT(backward->get_data()[i], DoubleNear(expected->get_data()[i], precision));
    }
  }
}

auto values_3 = make_shared<encap_t>(
  vector<double>{
    0.0,     third<double>(),     two_thirds<double>(),
    1.0, 1 + third<double>(), 1 + two_thirds<double>(),
    2.0, 2 + third<double>(), 2 + two_thirds<double>(),
    3.0, 3 + third<double>(), 3 + two_thirds<double>(),
    4.0, 4 + third<double>(), 4 + two_thirds<double>(),
    5.0, 5 + third<double>(), 5 + two_thirds<double>(),
    6.0, 6 + third<double>(), 6 + two_thirds<double>(),
    7.0, 7 + third<double>(), 7 + two_thirds<double>(),
    8.0, 8 + third<double>(), 8 + two_thirds<double>()
  }
);

auto values_4 = make_shared<encap_t>(
  vector<double>{
    0.0, 0.25, 0.5, 0.75,
    1.0, 1.25, 1.5, 1.75,
    2.0, 2.25, 2.5, 2.75,
    3.0, 3.25, 3.5, 3.75,
    4.0, 4.25, 4.5, 4.75,
    5.0, 5.25, 5.5, 5.75,
    6.0, 6.25, 6.5, 6.75,
    7.0, 7.25, 7.5, 7.75,
    8.0, 8.25, 8.5, 8.75,
    9.0, 9.25, 9.5, 9.75,
    10.0, 10.25, 10.5, 10.75,
    11.0, 11.25, 11.5, 11.75,
    12.0, 12.25, 12.5, 12.75,
    13.0, 13.25, 13.5, 13.75,
    14.0, 14.25, 14.5, 14.75,
    15.0, 15.25, 15.5, 15.75,
  }
);

const vector<double> generate_v8() {
  vector<double> v(512);
  for (size_t z = 0; z < 8; ++z) {
    for (size_t y = 0; y < 8; ++y) {
      for (size_t x = 0; x < 8; ++x) {
        v[pfasst::linearized_index(make_tuple(z, y, x), 8)] = z + y + (x/7);
      }
    }
  }
  return v;
}
auto values_8 = make_shared<encap_t>(generate_v8());

INSTANTIATE_TEST_CASE_P(DFT,
                        DiscreteFastFourierTransform,
                        ::testing::Values(
                            values_3
                          , values_4
                          , values_8
                        ));


TEST_MAIN()
