#include "fixtures/test_helpers.hpp"
using ::testing::DoubleNear;
using ::testing::NotNull;

#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <vector>
using namespace std;

#include <pfasst/contrib/fft.hpp>
using pfasst::contrib::FFT;

#include <pfasst/encap/vector.hpp>
using encap_t = pfasst::encap::VectorEncapsulation<double, double, 1>;
using fft_t = FFT<encap_t>;

#include <pfasst/logging.hpp>
#include <pfasst/globals.hpp>


using FFTTypes = ::testing::Types<fft_t>;
INSTANTIATE_TYPED_TEST_CASE_P(FFT1D, Concepts, FFTTypes);

class Interface
  : public ::testing::Test
{
  protected:
    fft_t fft;
};

TEST_F(Interface, query_z_pointer_for_specific_num_dofs)
{
  complex<double>* z_ptr = fft.get_workspace({{1}})->z;
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

    vector<complex<double>> eval_base_function(const shared_ptr<encap_t> vec, const size_t k)
    {
      const size_t ndofs = vec->get_data().size();
      vector<complex<double>> result(ndofs);

      transform(vec->get_data().cbegin(), vec->get_data().cend(),
                result.begin(),
                [ndofs, k](const double& t) {
                  return exp(complex<double>(0.0, 1.0) * (TWO_PI / double(ndofs) * k * t));
                });

      return result;
    }

    vector<double> two_pi_k_t(const shared_ptr<encap_t> vec, const size_t& k)
    {
      const size_t ndofs = vec->get_data().size();
      vector<double> result(ndofs);

      transform(vec->get_data().cbegin(), vec->get_data().cend(),
                result.begin(),
                [k](const double& t) {
                  return cos(TWO_PI * k * t);
                });

      return result;
    }
};

TEST_P(DiscreteFastFourierTransform, forward_transform)
{
  size_t ndofs = values->get_data().size();
  for (size_t k = 0; k < ndofs; ++k) {
    const double precision = k * ndofs * numeric_limits<double>::epsilon();

    vector<complex<double>> forward(ndofs);
    auto* fft_forward = fft.forward(make_shared<encap_t>(two_pi_k_t(values, k)));
    for (size_t i = 0; i < ndofs; ++i) {
      forward[i] = fft_forward[i];
    }

    for (size_t i = 0; i < ndofs; ++i) {
      EXPECT_THAT(forward[i].imag(), DoubleNear(0.0, precision));

      if (i != k && i != (ndofs - k)) {
        EXPECT_THAT(forward[i].real(), DoubleNear(0.0, precision));
      } else if (ndofs % 2 == 0 && (i == 0 || i == ndofs / 2)) {
          EXPECT_THAT(forward[i].real(), DoubleNear(ndofs, precision));
      } else {
        // TODO: test the other cases where entries are of ndofs/2
      }
    }
  }
}


TEST_P(DiscreteFastFourierTransform, backward_transform)
{
  size_t ndofs = values->get_data().size();
  for (size_t k = 0; k < ndofs; ++k) {
    const double precision = k * ndofs * numeric_limits<double>::epsilon();

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

auto values_3 = make_shared<encap_t>(vector<double>{0.0, (1.0/3.0), (2.0/3.0)});
auto values_4 = make_shared<encap_t>(vector<double>{0.0, 0.25, 0.5, 0.75});
auto values_5 = make_shared<encap_t>(vector<double>{0.0, 0.2, 0.4, 0.6, 0.8});
auto values_10 = make_shared<encap_t>(vector<double>{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9});

INSTANTIATE_TEST_CASE_P(DFT,
                        DiscreteFastFourierTransform,
                        ::testing::Values(values_4, values_5, values_10));


TEST_MAIN()
