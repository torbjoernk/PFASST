#ifndef _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
#define _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_

#include <memory>
#include <type_traits>
#include <vector>
using std::shared_ptr;
using std::vector;

#include <pfasst/sweeper/imex.hpp>
#include <pfasst/contrib/fft.hpp>


namespace pfasst
{
  namespace examples
  {
    /**
     * @defgroup Heat2D Heat2D
     * @ingroup Examples
     */

    /**
     * @ingroup Heat2D
     */
    namespace heat2d
    {
      /**
       * @ingroup Heat2D
       */
      template<
        class SweeperTrait,
        typename Enabled = void
      >
      class Heat2D
        : public IMEX<SweeperTrait, Enabled>
      {
        static_assert(std::is_same<
                        typename SweeperTrait::encap_t::traits::dim_t,
                        std::integral_constant<size_t, 2>
                      >::value,
                      "Heat2D Sweeper requires 2D data structures");

        public:
          using traits = SweeperTrait;

          static void init_opts();

        private:
          using spatial_t = typename traits::spatial_t;

          typename traits::time_t                        _t0{0.0};
          spatial_t                                      _nu{0.2};
          pfasst::contrib::FFT<typename traits::encap_t> _fft;
          vector<vector<spatial_t>>                      _lap;

        protected:
          virtual shared_ptr<typename SweeperTrait::encap_t>
          evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                            const shared_ptr<typename SweeperTrait::encap_t> u) override;

          virtual shared_ptr<typename SweeperTrait::encap_t>
          evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                            const shared_ptr<typename SweeperTrait::encap_t> u) override;

          virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                      shared_ptr<typename SweeperTrait::encap_t> u,
                                      const typename SweeperTrait::time_t& t,
                                      const typename SweeperTrait::time_t& dt,
                                      const shared_ptr<typename SweeperTrait::encap_t> rhs) override;

          virtual vector<shared_ptr<typename SweeperTrait::encap_t>>
          compute_error(const typename SweeperTrait::time_t& t);

          virtual vector<shared_ptr<typename SweeperTrait::encap_t>>
          compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                 const typename SweeperTrait::time_t& t);

        public:
          explicit Heat2D(const size_t ndofs);
          Heat2D(const Heat2D<SweeperTrait, Enabled>& other) = default;
          Heat2D(Heat2D<SweeperTrait, Enabled>&& other) = default;
          virtual ~Heat2D() = default;
          Heat2D<SweeperTrait, Enabled>& operator=(const Heat2D<SweeperTrait, Enabled>& other) = default;
          Heat2D<SweeperTrait, Enabled>& operator=(Heat2D<SweeperTrait, Enabled>&& other) = default;

          virtual void set_options() override;

          virtual shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t);

          virtual void post_step() override;

          virtual bool converged(const bool pre_check) override;
          virtual bool converged() override;

          size_t get_num_dofs() const;
      };
    }  // ::pfasst::examples::heat2d
  }  // ::pfasst::examples
}  // ::pfasst

#include "heat2d_sweeper_impl.hpp"

#endif  // _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
