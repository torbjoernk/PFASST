#ifndef _PFASST__ENCAP__TRAITS_HPP_
#define _PFASST__ENCAP__TRAITS_HPP_

#include <type_traits>
#include <vector>
using std::vector;

#include "pfasst/globals.hpp"


namespace pfasst
{
  namespace encap
  {
    struct encap_data_tag
    {};

    struct vector_encap_tag
      : public encap_data_tag
    {};

    struct eigen3_encap_tag
      : public encap_data_tag
    {};

    /**
     * Type Traits for encapsulation of user data types.
     *
     * @tparam TimePrecision    the time precision, e.g. precision of the integration nodes
     * @tparam SpatialPrecision the spatial data precision
     * @tparam DataT            the actual data type encapsulated
     *
     * @ingroup Traits
     */
    template<
      class TimePrecision,
      class SpatialPrecision,
      size_t Dim,
      class DataT,
      class... Ts
    >
    struct encap_traits
    {
      //! public member type for the time precision
      using time_t = TimePrecision;

      //! public member type for the spatial precision
      using spatial_t = SpatialPrecision;

      //! public member type for the encapsulated data type
      using data_t = DataT;

      using tag_t = encap_data_tag;

      using dim_t = std::integral_constant<size_t, Dim>;

      static constexpr size_t DIM = Dim;
    };


    /**
     * Spatialized Type Traits for encapsulation of std::vector.
     *
     * @tparam TimePrecision    the time precision, e.g. precision of the integration nodes
     * @tparam SpatialPrecision the spatial data precision
     *
     * @ingroup Traits
     */
    template<
      class TimePrecision,
      class SpatialPrecision,
      size_t Dim
    >
    struct vector_encap_traits
      : public encap_traits<TimePrecision, SpatialPrecision, Dim, vector<SpatialPrecision>>
    {
      using time_t = TimePrecision;
      using spatial_t = SpatialPrecision;
      using data_t = vector<spatial_t>;
      using tag_t = vector_encap_tag;
      using dim_t = std::integral_constant<size_t, Dim>;
      static constexpr size_t  DIM = Dim;
    };


    /**
     * Spatialized Type Traits for encapsulation of std::vector.
     *
     * @tparam TimePrecision    the time precision, e.g. precision of the integration nodes
     * @tparam SpatialPrecision the spatial data precision
     *
     * @ingroup Traits
     */
    template<
      class TimePrecision,
      class SpatialPrecision,
      size_t Dim
    >
    struct eigen3_encap_traits
      : public encap_traits<TimePrecision, SpatialPrecision, Dim, EigenVector<SpatialPrecision>>
    {
      using time_t = TimePrecision;
      using spatial_t = SpatialPrecision;
      using data_t = EigenVector<spatial_t>;
      using tag_t = eigen3_encap_tag;
      using dim_t = std::integral_constant<size_t, Dim>;
      static constexpr size_t  DIM = Dim;
    };
  }  // ::pfasst::encap
}  // ::pfasst

#endif  // _PFASST__ENCAP__TRAITS_HPP_
