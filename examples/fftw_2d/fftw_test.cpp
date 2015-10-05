#include <algorithm>
#include <cmath>
#include <cstdio>
#include <complex>
#include <iostream>
#include <memory>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <boost/math/constants/constants.hpp>
#include <leathers/pop>
using boost::math::constants::pi;
using boost::math::constants::two_pi;

#include "pfasst/util.hpp"
#include "pfasst/contrib/fft.hpp"
#include "pfasst/encap/vector.hpp"
using Vector = pfasst::encap::VectorEncapsulation<double, double, 2>;
using FFT = pfasst::contrib::FFT<Vector>;

double func_sin_cos(const double& x, const double& y)
{
  return sin(x) + cos(y);
}

double func_e_x(const double& x, const double& y)
{
  return exp(-1.0 * (pow(x - pi<double>(), 2) + pow(y - pi<double>(), 2)));
}

#define func func_sin_cos
// #define func func_e_x

void print_vec_2d(const shared_ptr<Vector>& vec, const pair<size_t, size_t>& dims) {
  for (size_t xi = 0; xi < dims.first; ++xi) {
    for (size_t yi = 0; yi < dims.second; ++yi) {
      printf("% 8.4f\t", vec->get_data()[xi * dims.first + yi]);
    }
    cout << endl;
  }
}

void print_arr_2d(const complex<double>* data, const pair<size_t, size_t>& dims) {
  for (size_t xi = 0; xi < dims.first; ++xi) {
    for (size_t yi = 0; yi < dims.second; ++yi) {
      printf("(% 8.4f +% 8.4fi)\t", real(data[xi * dims.first + yi]), imag(data[xi * dims.first + yi]));
    }
    cout << endl;
  }
}

double compute_diff(const shared_ptr<Vector>& a, const shared_ptr<Vector>& b)
{
  assert(a->get_data().size() == b->get_data().size());
  const size_t nx = a->get_data().size();

  auto diff = make_shared<Vector>(nx);
  transform(a->get_data().cbegin(), a->get_data().cend(),
            b->get_data().cbegin(),
            diff->data().begin(),
            [](const double& ai, const double& bi) { return abs(ai - bi); } );
  const auto norm = *max_element(diff->get_data().cbegin(), diff->get_data().cend());
  cout << "difference: (max-norm=" << norm << ")" << endl;
  print_vec_2d(diff, make_pair(sqrt(nx), sqrt(nx)));
  return norm;
}

void interpolate_data(FFT& fft, const shared_ptr<Vector>& coarse, shared_ptr<Vector> fine)
{
  const size_t coarse_ndofs = coarse->get_data().size();
  const size_t fine_ndofs = fine->get_data().size();
  assert(coarse_ndofs > 0);
  assert(fine_ndofs >= coarse_ndofs);

  if (fine_ndofs == coarse_ndofs) {
    // do a shortcut
    cout << "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT" << endl;
    fine->data() = coarse->get_data();

  } else {
    assert(coarse_ndofs % 2 == 0 && fine_ndofs % 2 == 0);
    const size_t coarse_dim_dofs = sqrt(coarse_ndofs);
    const size_t fine_dim_dofs   = sqrt(fine_ndofs);

    complex<double> *coarse_z = fft.forward(coarse);
    cout << "coarse: " << endl;
    print_arr_2d(coarse_z, make_pair(coarse_dim_dofs, coarse_dim_dofs));
    cout << endl;

    complex<double> *fine_z = fft.get_workspace(fine->get_dimwise_num_dofs())->z;

    for (size_t i = 0; i < fine_ndofs; i++) {
      fine_z[i] = 0.0;
    }

    // FFTW is not normalized
    double c = 1.0 / (double)coarse_ndofs;

    // a few utility functions
    auto is_in_first_half = [&](const size_t i) {
      return (i < coarse_dim_dofs / 2);
    };
    auto get_fine_dim_back_index = [&](const size_t ci) {
      return fine_dim_dofs - (fine_dim_dofs / 4) + ci - coarse_dim_dofs / 2;
    };
    auto get_fine_dim_index = [&](const size_t ci) {
      return is_in_first_half(ci) ? ci : get_fine_dim_back_index(ci);
    };
    auto get_coarse_index = [&](const size_t yi, const size_t xi) {
      return yi * coarse_dim_dofs + xi;
    };
    auto get_fine_index = [&](const size_t yi, const size_t xi) {
      return yi * fine_dim_dofs + xi;
    };

    for (size_t yi = 0; yi < coarse_dim_dofs; ++yi) {
      // y is second dim (i.e. columns)
      const size_t fine_yi = get_fine_dim_index(yi);

      for (size_t xi = 0; xi < coarse_dim_dofs; ++xi) {
        // x is first dim (i.e. rows)
        const size_t fine_xi = get_fine_dim_index(xi);

        const size_t coarse_index = get_coarse_index(yi, xi);
        if (coarse_index >= coarse_ndofs ) cout << "coarse index too large: " << coarse_index << " >= " << coarse_ndofs << endl;
        assert(coarse_index < coarse_ndofs);
        const size_t fine_index = get_fine_index(fine_yi, fine_xi);
        if (fine_index >= fine_ndofs) cout << "fine index too large: " << fine_index << " >= " << fine_ndofs << endl;
        assert(fine_index < fine_ndofs);

        fine_z[fine_index] = c * coarse_z[coarse_index];
      }
    }

    print_arr_2d(fine_z, make_pair(fine_dim_dofs, fine_dim_dofs));

    fft.backward(fine);
  }
}

void restrict_data(FFT& fft, const shared_ptr<Vector>& fine, shared_ptr<Vector> coarse)
{
  const size_t coarse_ndofs = coarse->get_data().size();
  const size_t fine_ndofs = fine->get_data().size();
  assert(coarse_ndofs > 0);
  assert(fine_ndofs >= coarse_ndofs);

  if (fine_ndofs == coarse_ndofs) {
    // do a shortcut
    coarse->data() = fine->get_data();

  } else {
    const size_t coarse_dim_dofs = sqrt(coarse_ndofs);
    const size_t fine_dim_dofs   = sqrt(fine_ndofs);
    const size_t factor = fine_dim_dofs / coarse_dim_dofs;
    cout << "coarsining by factor=" << factor << endl;

    for (size_t yi = 0; yi < coarse_dim_dofs; ++yi) {
      for (size_t xi = 0; xi < coarse_dim_dofs; ++xi) {
        const size_t coarse_index = yi * coarse_dim_dofs + xi;
        assert(coarse_index < coarse_ndofs);
        const size_t fine_index = factor * (yi * fine_dim_dofs + xi);
        assert(fine_index < fine_ndofs);
        coarse->data()[coarse_index] = fine->get_data()[fine_index];
      }
    }
  }
}


int main(int argc, char** argv) {
  assert(argc == 2 && atoi(argv[1]) > 0);
  using pfasst::almost_zero;

  const size_t dim_size = atoi(argv[1]);
  cout << "Using " << dim_size << " points in both dimensions" << endl;

  const double cdx = two_pi<double>() / double(dim_size);
  cout << "using cdx=cdy=" << cdx << endl;

  auto vec = make_shared<Vector>(dim_size * dim_size);

  // initialize data
  for (size_t xi = 0; xi < dim_size; ++xi) {
    const double x = 0 + cdx * xi;
    for (size_t yi = 0; yi < dim_size; ++yi) {
      const double y = 0 + cdx * yi;
      vec->data()[xi * dim_size + yi] = func(x, y);
    }
  }

  cout << endl << "initial data:" << endl;
  print_vec_2d(vec, make_pair(dim_size, dim_size));

  FFT fft;

  cout << endl << "FFT space:" << endl;
  auto* z_ptr = fft.forward(vec);
  print_arr_2d(z_ptr, make_pair(dim_size, dim_size));

  auto back = make_shared<Vector>(dim_size * dim_size);
  cout << endl << "back transform:" << endl;
  fft.backward(back);

  auto scaled_back = make_shared<Vector>(dim_size * dim_size);
  scaled_back->scaled_add(1 / double(dim_size * dim_size), back);
  print_vec_2d(scaled_back, make_pair(dim_size, dim_size));

  const size_t fine_dim_size = dim_size * 2;
  cout << endl << "transfer: nc=" << dim_size << "->" << fine_dim_size << "=fc" << endl;
  auto fine = make_shared<Vector>(fine_dim_size * fine_dim_size);
  interpolate_data(fft, vec, fine);

  cout << "interpolated:" << endl;
  print_vec_2d(fine, make_pair(fine_dim_size, fine_dim_size));

  const double fdx = two_pi<double>() / double(fine_dim_size);
  cout << "using fdx=fdy=" << cdx << endl;
  auto fine_fun = make_shared<Vector>(fine_dim_size * fine_dim_size);
  for (size_t xi = 0; xi < fine_dim_size; ++xi) {
    const double x = 0 + fdx * xi;
    for (size_t yi = 0; yi < fine_dim_size; ++yi) {
      const double y = 0 + fdx * yi;
      fine_fun->data()[xi * fine_dim_size + yi] = func(x, y);
    }
  }
  cout << "fine values:" << endl;
  print_vec_2d(fine_fun, make_pair(fine_dim_size, fine_dim_size));
  const auto norm = compute_diff(fine, fine_fun);

  assert(norm < 1e-15);

  auto coarse = make_shared<Vector>(dim_size * dim_size);
  restrict_data(fft, fine_fun, coarse);
  cout << "restricted fine:" << endl;
  print_vec_2d(coarse, make_pair(dim_size, dim_size));
  const auto norm_c = compute_diff(coarse, vec);
  assert(norm_c < 1e-15);
}
