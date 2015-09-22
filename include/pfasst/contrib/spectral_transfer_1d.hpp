#ifndef _PFASST__TRANSFER__SPECTRAL_TRANSFER_1D_HPP_
#define _PFASST__TRANSFER__SPECTRAL_TRANSFER_1D_HPP_

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
      class TransferTraits
    >
    class SpectralTransfer<TransferTraits, typename enable_if<
                 is_same<
                   typename TransferTraits::fine_encap_traits::dim_type,
                   integral_constant<size_t, 1>
                 >::value
               >::type>
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
        virtual void interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_type> coarse,
                                      shared_ptr<typename TransferTraits::fine_encap_type> fine);

        virtual void restrict_data(const shared_ptr<typename TransferTraits::fine_encap_type> fine,
                                   shared_ptr<typename TransferTraits::coarse_encap_type> coarse);
    };
  }  // ::pfasst::contrib
}  // ::pfasst

#include "pfasst/contrib/spectral_transfer_1d_impl.hpp"

#endif  // _PFASST__TRANSFER__SPECTRAL_TRANSFER_1D_HPP_
