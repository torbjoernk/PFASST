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
using boost::math::constants::two_pi;

#include "pfasst/contrib/fft_2d.hpp"
using FFT = pfasst::contrib::FFT2D<double>;
#include "pfasst/encap/vector.hpp"
using Vector = pfasst::encap::VectorEncapsulation<double, double>;

double func_sin_cos(const double& x, const double& y)
{
  return sin(x) + cos(y);
}

void print_vec_2d(const shared_ptr<Vector> vec, const pair<size_t, size_t>& dims) {
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


int main(int argc, char** argv) {
  assert(argc == 2 && atoi(argv[1]) > 0);

  const size_t dim_size = atoi(argv[1]);
  cout << "Using " << dim_size << " points in both dimensions" << endl;

  const double dx = two_pi<double>() / double(dim_size);
  cout << "using dx=dy=" << dx << endl;

  auto vec = make_shared<Vector>(dim_size * dim_size);

  // initialize data
  for (size_t xi = 0; xi < dim_size; ++xi) {
    const double x = 0 + dx * xi;
    for (size_t yi = 0; yi < dim_size; ++yi) {
      const double y = 0 + dx * yi;
      vec->data()[xi * dim_size + yi] = func_sin_cos(x, y);
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
}
