#ifndef _PFASST__EXAMPLES__ADVEC_DIFF__ADVEC_DIFF_SWEEPER_HPP_
#define _PFASST__EXAMPLES__ADVEC_DIFF__ADVEC_DIFF_SWEEPER_HPP_

#include <complex>
#include <memory>
#include <type_traits>
#include <vector>
using std::shared_ptr;
using std::vector;

#include <pfasst/sweeper/imex.hpp>
#include <pfasst/contrib/fft.hpp>

// I'd really like to have these as static const variable templates but this is only possible since C++14 ... :-(
//! Default diffusivity coefficient.
#define DEFAULT_DIFFUSIVITY 0.02
//! Default velocity.
#define DEFAULT_VELOCITY    1.0


namespace pfasst
{
  namespace examples
  {
    /**
     * @defgroup AdvectionDiffusion Advection-Diffusion
     * @brief Classic advection-diffusion equation in 1D solved with FFT in space.
     * @details This is an example demonstration the core features of the PFASST++ framework with
     *   the help of the classic advection-diffusion equation on a unit line with periodic boundary
     *   conditions.
     *
     * @ingroup Examples
     */

    /**
     * @ingroup AdvectionDiffusion
     */
    namespace advec_diff
    {
      /**
       * Advection-Diffusion Sweeper based an the IMEX SDC scheme.
       *
       * @f[
       *   \frac{\partial \vec{u}}{\partial t} = \nu \Delta \vec{u} - v \nabla \vec{u}
       * @f]
       * with @f$ \nu @f$ being the diffusivity coefficient and @f$ v @f$ the velocity of the
       * advection.
       *
       * To implement the implicit-explicit (IMEX) Spectral Deferred Corrections (SDC) scheme, the
       * evaluation of this equation is split into the explicit part (@f$ - v \nabla \vec{u} @f$)
       * and implicit part (@f$ \nu \Delta \vec{u} @f$).
       *
       * Both, the Laplace operator and the gradient are discretized using a FFT.
       *
       * @see `AdvecDiff::init_opts()` for additional command line parameters
       *
       * @ingroup AdvectionDiffusion
       */
      template<
        class SweeperTrait,
        typename Enabled = void
      >
      class AdvecDiff
        : public IMEX<SweeperTrait, Enabled>
      {
        static_assert(std::is_same<
                        typename SweeperTrait::encap_traits::dim_t,
                        std::integral_constant<size_t, 1>
                      >::value,
                      "Advection-Diffusion Sweeper requires 1D data structures");

        public:
          //! @copybrief IMEX::traits
          using traits = SweeperTrait;

          /**
           * Initialize additional command line options.
           */
          static void init_opts();

        private:
          //! @name Problem Parameters
          //! @{

          //! Initial time of the problem.
          typename traits::time_t    _t0;
          //! Diffusivity Coefficient.
          typename traits::spatial_t _nu;
          //! Velocity.
          typename traits::spatial_t _v;
          //! @}

          //! @name Utility Cache
          //! @{

          //! Instance of the FFT workspace.
          pfasst::contrib::FFT<typename traits::encap_t>   _fft;
          //! FFT coefficients for the gradient discretization in Fourier space.
          vector<std::complex<typename traits::spatial_t>> _ddx;
          //! FFT coefficients for the Laplace discretization in Fourier space.
          vector<std::complex<typename traits::spatial_t>> _lap;
          //! @}

        protected:
          //! @name Problem Equation Evaluation
          //! @{
          /**
           * @copybrief IMEX::evaluate_rhs_expl
           *
           * @param[in] t  time point to evaluate @f$ F_E(\vec{u},t) @f$ at
           * @param[in] u  spatial data to be passed to @f$ F_E(\vec{u},t) @f$
           * @returns spatial values of function evaluation
           */
          virtual shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                                               const shared_ptr<typename SweeperTrait::encap_t> u) override;

          /**
           * @copybrief IMEX::evaluate_rhs_impl
           *
           * @param[in] t  time point to evaluate @f$ F_I(\vec{u},t) @f$ at
           * @param[in] u  spatial data to be passed to @f$ F_I(\vec{u},t) @f$
           * @returns spatial values of function evaluation
           */
          virtual shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                                               const shared_ptr<typename SweeperTrait::encap_t> u) override;

          /**
           * @copybrief IMEX::implicit_solve
           *
           * The implicit SDC equation
           * @f[
           *   \vec{u}_{n+1} - \Delta_{t_{n+1}} F_I(\vec{u}_{n+1}, t_{n+1}) = \vec{u}_n
           *    + \Delta_{t_{n+1}} \left( F_E(\vec{u}_{n+1},t_{n+1}) - F_E(\vec{u}_n, t_n)
           *    - F_I(\vec{u}_n,t_n) + Q F_n \right)
           * @f]
           * is solved for @f$ \vec{u}_{n+1} @f$ and @f$ F_I(\vec{u}_{n+1}, t_{n+1}) @f$ in Fourier
           * space using FFT.
           *
           * @param[in,out] f    values of implicit function evaluation @f$ F_I(\vec{u}_{n+1}, t_{n+1}) @f$
           * @param[in,out] u    desired spatial data @f$ \vec{u}_{n+1} @f$
           * @param[in]     t    _target_ time point @f$ t_{n+1} @f$
           * @param[in]     dt   width of the current time interval
           * @param[in]     rhs  right hand side w.r.t. to SDC as defined by equation above
           */
          virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                      shared_ptr<typename SweeperTrait::encap_t> u,
                                      const typename SweeperTrait::time_t& t,
                                      const typename SweeperTrait::time_t& dt,
                                      const shared_ptr<typename SweeperTrait::encap_t> rhs) override;
          //! @}

          //! @name Analysis
          //! @{
          /**
           * Compute error at all nodes for current time step.
           *
           * The error is computed for each quadrature node by component wise difference compared to
           * the exact solution.
           *
           * @param[in] t time step's initial time point
           *
           * @returns vector of absolute errors at each time node
           */
          virtual vector<shared_ptr<typename SweeperTrait::encap_t>> compute_error(const typename SweeperTrait::time_t& t);
          /**
           * Compute errors at all nodes for current time step relative to norm of solution.
           *
           * Scales the absolute error at each time node to the norm of the spatial solution at the
           * respective time node:
           *
           * @f[
           *   \vec{e}_{rel} = \frac{\vec{e}_{abs}}{\| \vec{u} \|_0}
           * @f]
           *
           * @param[in] error  absolute errors at all time nodes
           * @param[in] t      @f$ t_0 @f$ of current time step
           *
           * @returns vector of errors scaled to the norm of their respective spatial data
           */
          virtual vector<shared_ptr<typename SweeperTrait::encap_t>> compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                                                                            const typename SweeperTrait::time_t& t);
          //! @}

        public:
          //! @name Constructors, Destructors & Common Operators
          //! @{
          explicit AdvecDiff(const size_t& ndofs, const typename SweeperTrait::spatial_t& nu = DEFAULT_DIFFUSIVITY,
                             const typename SweeperTrait::spatial_t& v = DEFAULT_VELOCITY);
          AdvecDiff(const AdvecDiff<SweeperTrait, Enabled>& other) = default;
          AdvecDiff(AdvecDiff<SweeperTrait, Enabled>&& other) = default;
          virtual ~AdvecDiff() = default;
          AdvecDiff<SweeperTrait, Enabled>& operator=(const AdvecDiff<SweeperTrait, Enabled>& other) = default;
          AdvecDiff<SweeperTrait, Enabled>& operator=(AdvecDiff<SweeperTrait, Enabled>&& other) = default;
          //! @}

          //! @name Configuration and Setup
          //! @{
          virtual void set_options() override;
          //! @}

          //! @name Problem Implementation
          //! @{
          /**
           * Computes the exact solution at given time point @p t of the problem equation.
           *
           * @param[in] t  time point of the desired exact solution
           */
          virtual shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t);
          //! @}

          //! @name Overwritten Functions
          //! @{
          /**
           * @copybrief IMEX::post_step()
           * @copydoc IMEX::post_step()
           */
          virtual void post_step() override;

          /**
           * @copybrief IMEX::converged(const bool)
           * @copydoc IMEX::converged(const bool)
           */
          virtual bool converged(const bool pre_check) override;
          /**
           * @copybrief IMEX::converged()
           * @copydoc IMEX:converged()
           */
          virtual bool converged() override;
          //! @}

          //! @name Utilities
          //! @{
          size_t get_num_dofs() const;
          //! @}
      };
    }  // ::pfasst::examples::advec_diff
  }  // ::pfasst::examples
}  // ::pfasst

#include "advec_diff_sweeper_impl.hpp"

#endif  // _PFASST__EXAMPLES__ADVEC_DIFF__ADVEC_DIFF_SWEEPER_HPP_
