#ifndef _PFASST__EXAMPLES__ADVEC_DIFF__ADVEC_DIFF_SWEEPER_HPP_
#define _PFASST__EXAMPLES__ADVEC_DIFF__ADVEC_DIFF_SWEEPER_HPP_

#include <memory>
#include <vector>
using namespace std;

#include <pfasst/sweeper/imex.hpp>
#include <pfasst/contrib/fft.hpp>

// I'd really like to have these as static const variable templates but this is only possible since C++14 ... :-(
#define DEFAULT_DIFFUSIVITY 0.02
#define DEFAULT_VELOCITY    1.0


namespace pfasst
{
  namespace examples
  {
    namespace advec_diff
    {
      template<
        class SweeperTrait,
        typename Enabled = void
      >
      class AdvecDiff
        : public IMEX<SweeperTrait, Enabled>
      {
        static_assert(is_same<
                        typename SweeperTrait::encap_traits::dim_t,
                        integral_constant<size_t, 1>
                      >::value,
                      "Advection-Diffusion Sweeper requires 1D data structures");

        public:
          using traits = SweeperTrait;

          static void init_opts();

        private:
          typename traits::time_t    _t0;
          typename traits::spatial_t _nu;
          typename traits::spatial_t _v;

          pfasst::contrib::FFT<typename traits::encap_t> _fft;
          vector<complex<typename traits::spatial_t>>    _ddx;
          vector<complex<typename traits::spatial_t>>    _lap;

        protected:
          virtual shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                                                  const shared_ptr<typename SweeperTrait::encap_t> u) override;
          virtual shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                                                  const shared_ptr<typename SweeperTrait::encap_t> u) override;

          virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                      shared_ptr<typename SweeperTrait::encap_t> u,
                                      const typename SweeperTrait::time_t& t,
                                      const typename SweeperTrait::time_t& dt,
                                      const shared_ptr<typename SweeperTrait::encap_t> rhs) override;

          virtual vector<shared_ptr<typename SweeperTrait::encap_t>> compute_error(const typename SweeperTrait::time_t& t);
          virtual vector<shared_ptr<typename SweeperTrait::encap_t>> compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                                                                            const typename SweeperTrait::time_t& t);

        public:
          explicit AdvecDiff(const size_t& ndofs, const typename SweeperTrait::spatial_t& nu = DEFAULT_DIFFUSIVITY,
                             const typename SweeperTrait::spatial_t& v = DEFAULT_VELOCITY);
          AdvecDiff(const AdvecDiff<SweeperTrait, Enabled>& other) = default;
          AdvecDiff(AdvecDiff<SweeperTrait, Enabled>&& other) = default;
          virtual ~AdvecDiff() = default;
          AdvecDiff<SweeperTrait, Enabled>& operator=(const AdvecDiff<SweeperTrait, Enabled>& other) = default;
          AdvecDiff<SweeperTrait, Enabled>& operator=(AdvecDiff<SweeperTrait, Enabled>&& other) = default;

          virtual void set_options() override;

          virtual shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t);

          virtual void post_step() override;

          virtual bool converged(const bool pre_check) override;
          virtual bool converged() override;

          size_t get_num_dofs() const;
      };
    }  // ::pfasst::examples::advec_diff
  }  // ::pfasst::examples
}  // ::pfasst

#include "advec_diff_sweeper_impl.hpp"

#endif  // _PFASST__EXAMPLES__ADVEC_DIFF__ADVEC_DIFF_SWEEPER_HPP_
