#ifndef _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
#define _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_

#include <memory>
#include <vector>
using namespace std;

#include <pfasst/sweeper/imex.hpp>
#include <pfasst/contrib/fft.hpp>


namespace pfasst
{
  namespace examples
  {
    namespace heat2d
    {
      template<
        class SweeperTrait,
        typename Enabled = void
      >
      class Heat2D
        : public IMEX<SweeperTrait, Enabled>
      {
        static_assert(is_same<
                        typename SweeperTrait::encap_type::traits::dim_type,
                        integral_constant<size_t, 2>
                      >::value,
                      "Heat2D Sweeper requires 2D data structures");

        public:
          typedef          SweeperTrait         traits;
          typedef typename traits::encap_type   encap_type;
          typedef typename traits::time_type    time_type;
          typedef typename traits::spatial_type spatial_type;

          static void init_opts();

        private:
          time_type    _t0;
          spatial_type _nu;

          pfasst::contrib::FFT<encap_type> _fft;
          vector<vector<spatial_type>>     _lap;

        protected:
          virtual shared_ptr<typename SweeperTrait::encap_type> evaluate_rhs_expl(const typename SweeperTrait::time_type& t,
                                                                                  const shared_ptr<typename SweeperTrait::encap_type> u) override;
          virtual shared_ptr<typename SweeperTrait::encap_type> evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                                                  const shared_ptr<typename SweeperTrait::encap_type> u) override;

          virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_type> f,
                                      shared_ptr<typename SweeperTrait::encap_type> u,
                                      const typename SweeperTrait::time_type& t,
                                      const typename SweeperTrait::time_type& dt,
                                      const shared_ptr<typename SweeperTrait::encap_type> rhs) override;

          virtual vector<shared_ptr<typename SweeperTrait::encap_type>> compute_error(const typename SweeperTrait::time_type& t);
          virtual vector<shared_ptr<typename SweeperTrait::encap_type>> compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_type>>& error, const typename SweeperTrait::time_type& t);

        public:
          explicit Heat2D(const size_t& ndofs, const typename SweeperTrait::spatial_type& nu = 0.02);
          Heat2D(const Heat2D<SweeperTrait, Enabled>& other) = default;
          Heat2D(Heat2D<SweeperTrait, Enabled>&& other) = default;
          virtual ~Heat2D() = default;
          Heat2D<SweeperTrait, Enabled>& operator=(const Heat2D<SweeperTrait, Enabled>& other) = default;
          Heat2D<SweeperTrait, Enabled>& operator=(Heat2D<SweeperTrait, Enabled>&& other) = default;

          virtual void set_options() override;

          virtual shared_ptr<typename SweeperTrait::encap_type> exact(const typename SweeperTrait::time_type& t);

          virtual void post_step() override;

          virtual bool converged(const bool& pre_check = false) override;

          size_t get_num_dofs() const;
      };
    }  // ::pfasst::examples::heat2d
  }  // ::pfasst::examples
}  // ::pfasst

#include "heat2d_sweeper_impl.hpp"

#endif  // _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
