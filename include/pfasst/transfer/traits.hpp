#ifndef _PFASST__TRANSFER__TRAITS_HPP_
#define _PFASST__TRANSFER__TRAITS_HPP_


namespace pfasst
{
  template<
    class CoarseSweeper,
    class FineSweeper,
    int NumLevels,
    class... Ts
  >
  struct transfer_traits
  {
    using coarse_sweeper_t = CoarseSweeper;
    using coarse_sweeper_traits = typename coarse_sweeper_t::traits;
    using coarse_encap_traits = typename coarse_sweeper_traits::encap_traits;
    using coarse_encap_t = typename coarse_sweeper_traits::encap_t;
    using coarse_time_t = typename coarse_sweeper_traits::time_t;
    using coarse_spatial_t = typename coarse_sweeper_traits::spatial_t;

    using fine_sweeper_t = FineSweeper;
    using fine_sweeper_traits = typename fine_sweeper_t::traits;
    using fine_encap_traits = typename fine_sweeper_traits::encap_traits;
    using fine_encap_t = typename fine_sweeper_traits::encap_t;
    using fine_time_t = typename fine_sweeper_traits::time_t;
    using fine_spatial_t = typename fine_sweeper_traits::spatial_t;

    static const size_t num_levels = NumLevels;
  };
}  // ::pfasst

#endif  // _PFASST__TRANSFER__TRAITS_HPP_
