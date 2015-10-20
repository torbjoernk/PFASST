#ifndef _PFASST__EXAMPLES__HEAD3D__HEAD3D_SPACE_PARALLEL_SWEEPER_HPP_
#define _PFASST__EXAMPLES__HEAD3D__HEAD3D_SPACE_PARALLEL_SWEEPER_HPP_

#include <memory>
#include <vector>
using namespace std;

#include <pfasst/sweeper/imex.hpp>
#include <pfasst/encap/traits.hpp>
#include <pfasst/contrib/fft.hpp>
#include <pfasst/comm/communicator.hpp>


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
                        typename SweeperTrait::encap_t::traits::dim_t,
                        integral_constant<size_t, 3>
                      >::value,
                      "Heat3D Sweeper requires 3D data structures");
        static_assert(is_same<
                        typename SweeperTrait::encap_t::traits::tag_t,
                        encap::vector_encap_tag
                      >::value,
                      "Heat3D Sweeper works only with std::vector Encapsulations");
        static_assert(is_base_of<
                        comm::Communicator<comm::spatial_communicator_tag>,
                        typename SweeperTrait::spatial_comm_t
                      >::value,
                      "space-parallel Heat3D requires spatial communicator");

        public:
          using traits = SweeperTrait;

          static void init_opts();

        private:
          using spatial_t = typename traits::spatial_t;

          typename traits::time_t                            _t0{0.0};
          spatial_t                                          _nu{0.2};
          pfasst::contrib::FFT<typename traits::encap_t>     _fft;
          vector<vector<vector<spatial_t>>>                  _lap;
          shared_ptr<typename traits::spatial_comm_t>        _space_comm;

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
          explicit Heat3D(const size_t ndofs);
          Heat3D(const Heat3D<SweeperTrait, Enabled>& other) = default;
          Heat3D(Heat3D<SweeperTrait, Enabled>&& other) = default;
          virtual ~Heat3D() = default;
          Heat3D<SweeperTrait, Enabled>& operator=(const Heat3D<SweeperTrait, Enabled>& other) = default;
          Heat3D<SweeperTrait, Enabled>& operator=(Heat3D<SweeperTrait, Enabled>&& other) = default;

          virtual void set_options() override;

          virtual       shared_ptr<typename SweeperTrait::spatial_comm_t>& spatial_communicator();
          virtual const shared_ptr<typename SweeperTrait::spatial_comm_t>& get_spatial_communicator() const;

          virtual shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t);

          virtual void post_step() override;

          virtual bool converged(const bool pre_check) override;
          virtual bool converged() override;

          size_t get_num_dofs() const;
      };
    }  // ::pfasst::examples::heat3d
  }  // ::pfasst::examples
}  // ::pfasst

#include "heat3d_space_parallel_sweeper_impl.hpp"

#endif  // _PFASST__EXAMPLES__HEAD3D__HEAD3D_SPACE_PARALLEL_SWEEPER_HPP_
