#ifndef _PFASST__TRANSFER__SPECTRAL_TRANSFER_HPP_
#define _PFASST__TRANSFER__SPECTRAL_TRANSFER_HPP_

#include "pfasst/transfer/polynomial.hpp"

#include <memory>
#include <vector>
using namespace std;

#include "pfasst/quadrature.hpp"
#include "pfasst/contrib/fft.hpp"


namespace pfasst
{
  namespace contrib
  {
    template<
      class TransferTraits,
      typename Enabled = void
    >
    class SpectralTransfer
      : public PolynomialTransfer<TransferTraits>
    {
      public:
        typedef TransferTraits traits;

        typedef typename traits::coarse_encap_traits coarse_encap_traits;
        typedef typename traits::coarse_encap_type coarse_encap_type;
        typedef typename traits::coarse_time_type coarse_time_type;
        typedef typename traits::coarse_spatial_type coarse_spatial_type;

        typedef typename traits::fine_encap_traits fine_encap_traits;
        typedef typename traits::fine_encap_type fine_encap_type;
        typedef typename traits::fine_time_type fine_time_type;
        typedef typename traits::fine_spatial_type fine_spatial_type;

      protected:
        pfasst::contrib::FFT<fine_encap_type> fft;

      public:
        SpectralTransfer() = default;
        SpectralTransfer(const SpectralTransfer<TransferTraits, Enabled> &other) = default;
        SpectralTransfer(SpectralTransfer<TransferTraits, Enabled> &&other) = default;
        virtual ~SpectralTransfer() = default;
        SpectralTransfer<TransferTraits, Enabled>& operator=(const SpectralTransfer<TransferTraits, Enabled> &other) = default;
        SpectralTransfer<TransferTraits, Enabled>& operator=(SpectralTransfer<TransferTraits, Enabled> &&other) = default;

        virtual void interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_type> coarse,
                                      shared_ptr<typename TransferTraits::fine_encap_type> fine);

        virtual void restrict_data(const shared_ptr<typename TransferTraits::fine_encap_type> fine,
                                   shared_ptr<typename TransferTraits::coarse_encap_type> coarse);
    };
  }  // ::pfasst::contrib
}  // ::pfasst

#include "pfasst/contrib/spectral_transfer_1d.hpp"
#include "pfasst/contrib/spectral_transfer_2d.hpp"

#endif  // _PFASST__TRANSFER__SPECTRAL_TRANSFER_HPP_
