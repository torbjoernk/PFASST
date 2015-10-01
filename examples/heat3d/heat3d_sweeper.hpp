#ifndef _PFASST__EXAMPLES__HEAD3D__HEAD3D_SWEEPER_HPP_
#define _PFASST__EXAMPLES__HEAD3D__HEAD3D_SWEEPER_HPP_

#include <memory>
#include <vector>
using namespace std;

#include <pfasst/sweeper/imex.hpp>
#include <pfasst/encap/traits.hpp>
#include <pfasst/contrib/fft.hpp>


namespace pfasst
{
  namespace examples
  {
    namespace heat3d
    {
      template<
        class SweeperTrait,
        typename Enabled = void
      >
      class Heat3D
        : public IMEX<SweeperTrait, Enabled>
      {
        static_assert(is_same<
                        typename SweeperTrait::encap_type::traits::dim_type,
                        integral_constant<size_t, 3>
                      >::value,
                      "Heat3D Sweeper requires 3D data structures");
        static_assert(is_same<
                        typename SweeperTrait::encap_type::traits::tag_type,
                        encap::vector_encap_tag
                      >::value,
                      "Heat3D Sweeper works only with std::vector Encapsulations");

        public:
          typedef          SweeperTrait         traits;
          typedef typename traits::encap_type   encap_type;
          typedef typename traits::time_type    time_type;
          typedef typename traits::spatial_type spatial_type;

          static void init_opts();

        private:
          time_type    _t0;
          spatial_type _nu;

          pfasst::contrib::FFT<encap_type>     _fft;
          vector<vector<vector<spatial_type>>> _lap;

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
          explicit Heat3D(const size_t& ndofs, const typename SweeperTrait::spatial_type& nu = 0.02);
          Heat3D(const Heat3D<SweeperTrait, Enabled>& other) = default;
          Heat3D(Heat3D<SweeperTrait, Enabled>&& other) = default;
          virtual ~Heat3D() = default;
          Heat3D<SweeperTrait, Enabled>& operator=(const Heat3D<SweeperTrait, Enabled>& other) = default;
          Heat3D<SweeperTrait, Enabled>& operator=(Heat3D<SweeperTrait, Enabled>&& other) = default;

          virtual void set_options() override;

          virtual shared_ptr<typename SweeperTrait::encap_type> exact(const typename SweeperTrait::time_type& t);

          virtual void post_step() override;

          virtual bool converged(const bool& pre_check = false) override;

          size_t get_num_dofs() const;
      };
    }  // ::pfasst::examples::heat3d
  }  // ::pfasst::examples
}  // ::pfasst

#include "heat3d_sweeper_impl.hpp"

#endif  // _PFASST__EXAMPLES__HEAD3D__HEAD3D_SWEEPER_HPP_
