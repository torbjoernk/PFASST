#include "pfasst/contrib/spectral_transfer.hpp"

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
    /**
     * @internal
     * The interpolation is done in the Fourier space based on the data representation of FFTW3.
     *
     * The positive frequencies are located in the top-left corner of the 2D matrix, while the
     * negative frequencies are in the bottom-right corner:
     *
     * @verbatim
     * + + . .
     * + + . .
     * . . - -
     * . . - -
     * @endverbatim
     *
     * Interpolation is then simply inserting a "plus" of zeros in the center of the matrix of
     * frequencies keeping the general order of the frequencies:
     *
     * @verbatim
     * + + 0 0 0 0 . .
     * + + 0 0 0 0 . .
     * 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0 0
     * . . 0 0 0 0 - -
     * . . 0 0 0 0 - -
     * @endverbatim
     * @endinterl
     */
    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits,
      typename enable_if<
                 is_same<
                   typename TransferTraits::fine_encap_traits::dim_type,
                   integral_constant<size_t, 2>
                 >::value
               >::type>::interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_type> coarse,
                                          shared_ptr<typename TransferTraits::fine_encap_type> fine)
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
        complex<fine_spatial_type> *coarse_z = this->fft.forward(coarse);
        complex<fine_spatial_type> *fine_z = this->fft.get_workspace(fine->get_dimwise_num_dofs())->z;

        const size_t coarse_dim_dofs = sqrt(coarse_ndofs);
        const size_t fine_dim_dofs   = sqrt(fine_ndofs);

        if (fine_dim_dofs != coarse_dim_dofs * 2) {
          ML_CLOG(FATAL, "TRANS", "FFTW based interpolation in 2D only for coarsening factor of 2");
          throw runtime_error("unsupported coarsening factor for FFTW interpolation");
        }

        for (size_t i = 0; i < fine_ndofs; i++) {
          fine_z[i] = 0.0;
        }

        // FFTW is not normalized
        double c = 1.0 / (double)coarse_ndofs;

        for (size_t yi = 0; yi < coarse_dim_dofs; ++yi) {
          // y is second dim (i.e. columns)
          for (size_t xi = 0; xi < coarse_dim_dofs; ++xi) {
            // x is first dim (i.e. rows)
            const size_t coarse_index = yi * coarse_dim_dofs + xi;
            assert(coarse_index < coarse_ndofs);

            if (yi < coarse_dim_dofs / 2 && xi < coarse_dim_dofs / 2) {
              // positive frequencies (in top-left corner)
              const size_t fine_index = yi * fine_dim_dofs + xi;
              assert(fine_index < fine_ndofs);
              fine_z[fine_index] = c * coarse_z[coarse_index];

            } else if (yi < coarse_dim_dofs / 2 && xi >= coarse_dim_dofs / 2) {
              // x-negative, y-positive frequencies (in top-right corner)
              const size_t fine_tail_col = fine_dim_dofs - (fine_dim_dofs / 4) + xi - coarse_dim_dofs / 2;
              const size_t fine_index = yi * fine_dim_dofs + fine_tail_col;
              assert(fine_index < fine_ndofs);
              fine_z[fine_index] = c * coarse_z[coarse_index];

            } else if (yi >= coarse_dim_dofs / 2 && xi < coarse_dim_dofs / 2) {
              // x-positive, y-negative frequencies (in bottom-left corner)
              const size_t fine_tail_row = fine_dim_dofs - (fine_dim_dofs / 4) + yi - coarse_dim_dofs / 2;
              const size_t fine_index = fine_tail_row * fine_dim_dofs + xi;
              assert(fine_index < fine_ndofs);
              fine_z[fine_index] = c * coarse_z[coarse_index];

            } else if (yi >= coarse_dim_dofs / 2 && xi >= coarse_dim_dofs / 2) {
              // negative frequencies (bottom-right corner)
              const size_t fine_tail_row = fine_dim_dofs - (fine_dim_dofs / 4) + yi - coarse_dim_dofs / 2;
              const size_t fine_tail_col = fine_dim_dofs - (fine_dim_dofs / 4) + xi - coarse_dim_dofs / 2;
              const size_t fine_index = fine_tail_row * fine_dim_dofs + fine_tail_col;
              assert(fine_index < fine_ndofs);
              fine_z[fine_index] = c * coarse_z[coarse_index];

            } else {
              // fine center null-plus
              continue;
            }
          }
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
                   typename TransferTraits::fine_encap_traits::dim_type,
                   integral_constant<size_t, 2>
                 >::value
               >::type>::restrict_data(const shared_ptr<typename TransferTraits::fine_encap_type> fine,
                                       shared_ptr<typename TransferTraits::coarse_encap_type> coarse)
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
        const size_t coarse_dim_dofs = sqrt(coarse_ndofs);
        const size_t fine_dim_dofs   = sqrt(fine_ndofs);
        const size_t factor = fine_dim_dofs / coarse_dim_dofs;

        if (fine_dim_dofs != coarse_dim_dofs * 2) {
          ML_CLOG(FATAL, "TRANS", "FFTW based interpolation in 2D only for coarsening factor of 2");
          throw runtime_error("unsupported coarsening factor for FFTW interpolation");
        }

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
  }  // ::pfasst::contrib
}  // ::pfasst
