#ifndef _PFASST__SWEEPER__IMEX_HPP_
#define _PFASST__SWEEPER__IMEX_HPP_

#include <memory>
#include <vector>
using std::shared_ptr;
using std::vector;

#include "pfasst/sweeper/sweeper.hpp"


namespace pfasst
{
  template<
    class SweeperTrait,
    typename Enabled = void
  >
  class IMEX
    : public Sweeper<SweeperTrait, Enabled>
  {
    public:
      using traits = SweeperTrait;

    protected:
      Matrix<typename traits::time_t> _q_delta_expl;
      Matrix<typename traits::time_t> _q_delta_impl;

      //! size = #nodes + 1
      vector<shared_ptr<typename traits::encap_t>> _q_integrals;
      vector<shared_ptr<typename traits::encap_t>> _expl_rhs;
      vector<shared_ptr<typename traits::encap_t>> _impl_rhs;

      size_t _num_expl_f_evals;
      size_t _num_impl_f_evals;
      size_t _num_impl_solves;

      virtual void integrate_end_state(const typename SweeperTrait::time_t& dt) override;
      virtual void compute_residuals(const bool& only_last) override;
      virtual void compute_residuals() override;
      virtual void initialize() override;

      virtual shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                                              const shared_ptr<typename SweeperTrait::encap_t> u);
      virtual shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                                              const shared_ptr<typename SweeperTrait::encap_t> u);

      virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                  shared_ptr<typename SweeperTrait::encap_t> u,
                                  const typename SweeperTrait::time_t& t,
                                  const typename SweeperTrait::time_t& dt,
                                  const shared_ptr<typename SweeperTrait::encap_t> rhs);

      virtual void compute_delta_matrices();

    public:
      explicit IMEX();
      IMEX(const IMEX<SweeperTrait, Enabled>& other) = default;
      IMEX(IMEX<SweeperTrait, Enabled>&& other) = default;
      virtual ~IMEX() = default;
      IMEX<SweeperTrait, Enabled>& operator=(const IMEX<SweeperTrait, Enabled>& other) = default;
      IMEX<SweeperTrait, Enabled>& operator=(IMEX<SweeperTrait, Enabled>&& other) = default;

      virtual void setup() override;

      virtual void pre_predict() override;
      /**
       * first order Euler method
       */
      virtual void predict() override;
      virtual void post_predict() override;

      virtual void pre_sweep() override;
      virtual void sweep() override;
      virtual void post_sweep() override;

      virtual void post_step() override;
      virtual void advance(const size_t& num_steps) override;
      virtual void advance() override;
      virtual void reevaluate(const bool initial_only) override;
      virtual void reevaluate() override;
      virtual vector<shared_ptr<typename SweeperTrait::encap_t>> integrate(const typename SweeperTrait::time_t& dt) override;
  };
}  // ::pfasst

#include "pfasst/sweeper/imex_impl.hpp"

#endif  // _PFASST__SWEEPER__IMEX_HPP_
