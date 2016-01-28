/**
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/fft_manager.hpp
 * @since v0.6.0
 */
#ifndef _EXAMPLES__ADVEC_DIFF__FFT_MANAGER_HPP_
#define _EXAMPLES__ADVEC_DIFF__FFT_MANAGER_HPP_

#include <map>
#include <memory>
using std::map;
using std::shared_ptr;


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      /**
       * Container to hold and query for FFTW workspaces.
       *
       * The @ref FFTManager holds all instances of @ref FFTWorkspace once queried through
       * FFTManager::get_workspace().
       *
       * @tparam WorkspaceT type of the @ref FFTWorkspace the manager instance is for
       * @ingroup AdvectionDiffusion
       */
      template<class WorkspaceT>
      class FFTManager
      {
        public:
          //! @{
          using workspace_t = WorkspaceT;
          //! @}

        protected:
          //! @{
          //! Storage of workspaces mapped to their number of DOFs
          map<size_t, shared_ptr<workspace_t>> _workspaces;
          //! @}

        public:
          //! @{
          FFTManager() = default;
          FFTManager(const FFTManager&) = default;
          FFTManager(FFTManager&&) = default;
          virtual ~FFTManager() = default;
          //! @}

          //! @{
          FFTManager& operator=(const FFTManager&) = default;
          FFTManager& operator=(FFTManager&&) = default;
          //! @}

          //! @{
          /**
           * Get the one @ref FFTWorkspace for given number of DOFs
           *
           * @param[in] ndofs number of degrees of freedom of FFTW workspace
           * @return shared pointer to the only workspace with given number of DOFs
           */
          shared_ptr<WorkspaceT> get_workspace(const size_t ndofs);
          //! @}

          //! @{
          /**
           * Finalizes cleanup of FFT resources
           *
           * Calls `WorkspaceT::finalize_cleanup()`.
           */
          static void finalize_cleanup();
          //! @}
      };
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#include "fft_manager_impl.hpp"

#endif // _EXAMPLES__ADVEC_DIFF__FFT_MANAGER_HPP_
