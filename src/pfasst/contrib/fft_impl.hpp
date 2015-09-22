#include "pfasst/contrib/fft.hpp"

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  namespace contrib
  {
    template<class Encapsulation>
    FFT<Encapsulation>::~FFT()
    {
      for (auto& x : workspaces) {
        shared_ptr<workspace> wk = std::get<1>(x);
        fftw_free(wk->wk);
        fftw_destroy_plan(wk->ffft);
        fftw_destroy_plan(wk->ifft);
      }
      workspaces.clear();
    }

    template<class Encapsulation>
    complex<typename Encapsulation::spatial_type>*
    FFT<Encapsulation>::forward(const shared_ptr<Encapsulation> x)
    {
      const auto ndofs = x->get_dimwise_num_dofs();
      auto wk = get_workspace(ndofs);
      for (size_t i = 0; i < x->get_total_num_dofs(); i++) {
        wk->z[i] = x->get_data()[i];
      }
      fftw_execute_dft(wk->ffft, wk->wk, wk->wk);
      return wk->z;
    }

    template<class Encapsulation>
    void
    FFT<Encapsulation>::backward(shared_ptr<Encapsulation> x)
    {
      const auto ndofs = x->get_dimwise_num_dofs();
      auto wk = get_workspace(ndofs);
      fftw_execute_dft(wk->ifft, wk->wk, wk->wk);
      for (size_t i = 0; i < x->get_total_num_dofs(); i++) {
        x->data()[i] = real(wk->z[i]);
      }
    }

    template<class Encapsulation>
    shared_ptr<typename FFT<Encapsulation>::workspace>
    FFT<Encapsulation>::get_workspace(const array<int, Encapsulation::traits::DIM>& ndofs)
    {
      const size_t DIM = Encapsulation::traits::DIM;
      const size_t total_ndofs = accumulate(ndofs.cbegin(), ndofs.cend(), 1, std::multiplies<int>());

      if (workspaces.find(ndofs) == workspaces.end()) {
        auto wk = make_shared<workspace>();
        wk->wk = fftw_alloc_complex(total_ndofs);
        wk->ffft = fftw_plan_dft(DIM, ndofs.data(), wk->wk, wk->wk, FFTW_FORWARD, FFTW_ESTIMATE);
        wk->ifft = fftw_plan_dft(DIM, ndofs.data(), wk->wk, wk->wk, FFTW_BACKWARD, FFTW_ESTIMATE);
        wk->z = reinterpret_cast<complex<typename Encapsulation::spatial_type>*>(wk->wk);
        workspaces.insert(pair<array<int, Encapsulation::traits::DIM>, shared_ptr<workspace>>(ndofs, wk));
      }

      return workspaces[ndofs];
    }
  }  // ::pfast::contrib
}  // ::pfasst
