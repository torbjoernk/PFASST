#ifndef _PFASST__SWEEPER__TRAITS_HPP_
#define _PFASST__SWEEPER__TRAITS_HPP_

#include "pfasst/encap/encapsulation.hpp"


namespace pfasst
{
  template<
    class EncapsulationTraits,
    class... Ts
  >
  struct sweeper_traits
  {
    using encap_traits = EncapsulationTraits;
    using encap_t = encap::Encapsulation<EncapsulationTraits>;
    using time_t = typename encap_traits::time_t;
    using spatial_t = typename encap_traits::spatial_t;
  };

  template<
    class EncapsulationTraits,
    class SpaceCommT,
    class... Ts
  >
  struct space_parallel_sweeper_traits
    : public sweeper_traits<EncapsulationTraits, SpaceCommT, Ts...>
  {
    using encap_traits = EncapsulationTraits;
    using encap_t = encap::Encapsulation<EncapsulationTraits>;
    using time_t = typename encap_traits::time_t;
    using spatial_t = typename encap_traits::spatial_t;
    using spatial_comm_t = SpaceCommT;
  };
}  // ::pfasst

#endif  // _PFASST__SWEEPER__TRAITS_HPP_
