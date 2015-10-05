#include "pfasst/contrib/spectral_transfer_1d.hpp"

#include <cassert>
#include <memory>
#include <vector>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/quadrature.hpp"


namespace pfasst
{
  namespace contrib
  {
    template<class TransferTraits>
    size_t
    SpectralTransfer<
      TransferTraits,
      typename enable_if<
                 is_same<
                   typename TransferTraits::fine_encap_traits::dim_t,
                   integral_constant<size_t, 1>
                 >::value
               >::type>::translate_index(const index_t& index, const index_t& extends) const
    {
      UNUSED(extends);
      return get<0>(index);
    }

    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits,
      typename enable_if<
                 is_same<
                   typename TransferTraits::fine_encap_traits::dim_t,
                   integral_constant<size_t, 1>
                 >::value
               >::type>::interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_t> coarse,
                                          shared_ptr<typename TransferTraits::fine_encap_t> fine)
    {
      ML_CVLOG(1, "TRANS", "interpolate data");

      const size_t coarse_ndofs = coarse->get_data().size();
      const size_t fine_ndofs = fine->get_data().size();
      assert(coarse_ndofs > 0);
      assert(fine_ndofs >= coarse_ndofs);

      if (fine_ndofs == coarse_ndofs) {
        // do a shortcut
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        fine->data() = coarse->get_data();

      } else {
        complex<typename traits::fine_spatial_t> *coarse_z = this->fft.forward(coarse);
        complex<typename traits::fine_spatial_t> *fine_z = this->fft.get_workspace(fine->get_dimwise_num_dofs())->z;

        for (size_t i = 0; i < fine_ndofs; i++) {
          fine_z[i] = 0.0;
        }

        // FFTW is not normalized
        typename traits::coarse_spatial_t c = 1.0 / coarse_ndofs;

        // positive frequencies
        for (size_t i = 0; i < coarse_ndofs / 2; i++) {
          fine_z[i] = c * coarse_z[i];
        }

        // negative frequencies (in backward order)
        for (size_t i = 1; i < coarse_ndofs / 2; i++) {
          fine_z[fine_ndofs - coarse_ndofs / 2 + i] = c * coarse_z[coarse_ndofs / 2 + i];
        }

        this->fft.backward(fine);
      }
    }

    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits,
      typename enable_if<
                 is_same<
                   typename TransferTraits::fine_encap_traits::dim_t,
                   integral_constant<size_t, 1>
                 >::value
               >::type>::restrict_data(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                       shared_ptr<typename TransferTraits::coarse_encap_t> coarse)
    {
      ML_CVLOG(1, "TRANS", "restrict data");

      const size_t coarse_ndofs = coarse->get_data().size();
      const size_t fine_ndofs = fine->get_data().size();
      assert(coarse_ndofs > 0);
      assert(fine_ndofs >= coarse_ndofs);

      if (fine_ndofs == coarse_ndofs) {
        // do a shortcut
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        coarse->data() = fine->get_data();

      } else {
        const size_t factor = fine_ndofs / coarse_ndofs;

        for (size_t i = 0; i < coarse_ndofs; i++) {
          coarse->data()[i] = fine->get_data()[factor * i];
        }
      }
    }
  }  // ::pfasst::contrib
}  // ::pfasst
