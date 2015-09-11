#include "pfasst/contrib/fft_2d.hpp"

namespace pfasst
{
  namespace contrib
  {
    template<class precision>
    FFT2D<precision>::~FFT2D()
    {
      for (auto& x : workspaces) {
        shared_ptr<workspace> wk = std::get<1>(x);
        fftw_free(wk->wk);
        fftw_destroy_plan(wk->ffft);
        fftw_destroy_plan(wk->ifft);
      }
      workspaces.clear();
    }

    template<class precision>
    shared_ptr<typename FFT2D<precision>::workspace>
    FFT2D<precision>::get_workspace(size_t ndofs)
    {
      if (workspaces.find(ndofs) == workspaces.end()) {
        auto wk = make_shared<workspace>();
        const size_t ndofs_dim = sqrt(ndofs);
        wk->wk = fftw_alloc_complex(ndofs);
        wk->ffft = fftw_plan_dft_2d(ndofs_dim, ndofs_dim, wk->wk, wk->wk, FFTW_FORWARD, FFTW_ESTIMATE);
        wk->ifft = fftw_plan_dft_2d(ndofs_dim, ndofs_dim, wk->wk, wk->wk, FFTW_BACKWARD, FFTW_ESTIMATE);
        wk->z = reinterpret_cast<complex<precision>*>(wk->wk);
        workspaces.insert(pair<size_t, shared_ptr<workspace>>(ndofs, wk));
      }

      return workspaces[ndofs];
    }

    template<class precision>
    template<class time_precision>
    complex<precision>*
    FFT2D<precision>::forward(const shared_ptr<encap::VectorEncapsulation<time_precision, precision>> x)
    {
      const size_t ndofs = x->get_data().size();
      auto wk = get_workspace(ndofs);
      for (size_t i = 0; i < ndofs; i++) {
        wk->z[i] = x->get_data()[i];
      }
      fftw_execute_dft(wk->ffft, wk->wk, wk->wk);
      return wk->z;
    }

    template<class precision>
    void
    FFT2D<precision>::backward(shared_ptr<encap::VectorEncapsulation<time_precision, precision>> x)
    {
      const size_t ndofs = x->get_data().size();
      auto wk = get_workspace(ndofs);
      fftw_execute_dft(wk->ifft, wk->wk, wk->wk);
      for (size_t i = 0; i < ndofs; i++) {
        x->data()[i] = real(wk->z[i]);
      }
    }
  }  // ::pfast::contrib
}  // ::pfasst
