/**
 * @ingroup Concepts
 * @file doc/source/concepts/fft_workspace.hpp
 * @since v0.6.0
 */
#include <complex>
using std::complex;


namespace pfasst
{
  /**
   * Specifications of Concepts used in PFASST
   *
   * @note This contains only the specifications of concepts and are not available as code.
   *
   * @ingroup Concepts
   */
  namespace concepts
  {
    /**
     * Concept of a FFT workspace
     *
     * A FFT workspace, manageable by @ref FFTManager, provides a wrapper around calls to FFT
     * and persists variables and memory between calls to FFT.
     *
     * A @ref FFTWorkspace is expected to conform to the RAII principle.
     *
     * Usually, the initial setup of FFT for a given number of degrees of freedom is done on
     * construction of a @ref FFTWorkspace, e.g. allocating memory for transformed values and
     * creating _plans_.
     *
     * Cleanup and freeing of this setup should be done on destruction.
     *
     * @tparam DataT Encapsulation type, e.g. pfasst::encap::VectorEncapsulation<double>;
     *               must provide public member type `DataT::value_type` for the type of a single
     *               DOF
     *
     * @note This is only the specification of a concept and not available as code.
     */
    template<class DataT>
    class FFTWorkspace
    {
      public:
        /**
         * Transforms given Encapsulation to Fourier space
         *
         * @param[in] x data in problem space to be transformed into Fourier space
         * @return pointer to transformed values
         */
        complex<typename DataT::value_type>* forward(const DataT& x);

        /**
         * Transformes data in Fourier space back to problem space
         *
         * Data in Fourier space as stored in `z_ptr()` are back-transformed into problem space
         * and stored in given `DataT` object, overwriting existing data.
         *
         * @param[in,out] x problem space data
         */
        void backward(DataT& x);

        /**
         * Number of degrees of Freedom of this workspace
         *
         * @return number of degrees of freedom
         */
        size_t size() const;

        /**
         * Accessor to data in Fourier space
         *
         * @return pointer to data in Fourier space
         */
        complex<DataT::value_type>* z_ptr();

        /**
         * Caller to FFT library specific cleanup routines
         *
         * Some FFT implementations require to call a cleanup function before program exit, e.g.
         * `fftw_cleanup()`.
         */
        static void finalize_cleanup();
    };
  }  // ::pfasst::concepts
}  // ::pfasst
