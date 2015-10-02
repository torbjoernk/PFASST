#include "fixtures/test_helpers.hpp"

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
using encap_t = pfasst::encap::VectorEncapsulation<double, double, 2>;
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
  complex<double>* z_ptr = fft.get_workspace({1, 1})->z;
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
      const size_t ndofs = vec->get_data().size();
      const size_t dim_ndofs = sqrt(ndofs);
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
  const size_t dim_ndofs = sqrt(ndofs);

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
    2.0, 2 + third<double>(), 2 + two_thirds<double>()
  }
);
auto values_4 = make_shared<encap_t>(
  vector<double>{
    0.0, 0.25, 0.5, 0.75,
    1.0, 1.25, 1.5, 1.75,
    2.0, 2.25, 2.5, 2.75,
    3.0, 3.25, 3.5, 3.75
  }
);
auto values_5 = make_shared<encap_t>(
  vector<double>{
    0.0, 0.2, 0.4, 0.6, 0.8,
    1.0, 1.2, 1.4, 1.6, 1.8,
    2.0, 2.2, 2.4, 2.6, 2.8,
    3.0, 3.2, 3.4, 3.6, 3.8,
    4.0, 4.2, 4.4, 4.6, 4.8
  }
);
auto values_10 = make_shared<encap_t>(
  vector<double>{
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
    2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
    3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
    4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9,
    5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9,
    6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9,
    7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9,
    8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9,
    9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9
  }
);

INSTANTIATE_TEST_CASE_P(DFT,
                        DiscreteFastFourierTransform,
                        ::testing::Values(values_3, values_4, values_5, values_10));


TEST_MAIN()
