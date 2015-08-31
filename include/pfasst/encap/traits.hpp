#ifndef _PFASST__ENCAP__TRAITS_HPP_
#define _PFASST__ENCAP__TRAITS_HPP_

#include <vector>
using namespace std;


namespace pfasst
{
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
    class SpatialPrecision
  >
  struct vector_encap_traits
    : public encap_traits<TimePrecision, SpatialPrecision, vector<SpatialPrecision>>
  {
    typedef TimePrecision        time_type;
    typedef SpatialPrecision     spatial_type;
    typedef vector<spatial_type> data_type;
  };

}  // ::pfasst

#endif  // _PFASST__ENCAP__TRAITS_HPP_
