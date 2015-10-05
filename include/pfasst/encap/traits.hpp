#ifndef _PFASST__ENCAP__TRAITS_HPP_
#define _PFASST__ENCAP__TRAITS_HPP_

#include <vector>
using namespace std;

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
  }  // ::pfasst::encap

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
    typedef TimePrecision    time_type;

    //! public member type for the spatial precision
    typedef SpatialPrecision spatial_type;

    //! public member type for the encapsulated data type
    typedef DataT            data_type;

    typedef integral_constant<size_t, Dim> dim_type;
    static constexpr size_t  DIM = Dim;

    typedef encap::encap_data_tag tag_type;
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
    typedef TimePrecision        time_type;
    typedef SpatialPrecision     spatial_type;
    typedef vector<spatial_type> data_type;
    typedef integral_constant<size_t, Dim> dim_type;
    static constexpr size_t  DIM = Dim;
    typedef encap::vector_encap_tag tag_type;
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
    size_t Dimension
  >
  struct eigen3_encap_traits
    : public encap_traits<TimePrecision, SpatialPrecision, Dimension, EigenVector<SpatialPrecision>>
  {
    typedef TimePrecision             time_type;
    typedef SpatialPrecision          spatial_type;
    typedef EigenVector<spatial_type> data_type;

    typedef integral_constant<size_t, Dimension> dim_type;
    static constexpr size_t  DIM = Dimension;

    typedef encap::eigen3_encap_tag tag_type;
  };
}  // ::pfasst

#endif  // _PFASST__ENCAP__TRAITS_HPP_
