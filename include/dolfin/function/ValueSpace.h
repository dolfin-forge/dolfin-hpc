
#ifndef __DOLFIN_FUNCTION_VALUE_SPACE_H_
#define __DOLFIN_FUNCTION_VALUE_SPACE_H_

#include <dolfin/common/types.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

template < size_t I = 1, size_t J = 1 >
class ValueSpace
{

public:
  ValueSpace()
  {
  }

  /// Return the rank of the value space
  static inline auto rank() -> size_t
  {
    return RANK;
  }

  /// Return the dimension of the value space for axis i
  static inline auto dim( size_t i ) -> size_t
  {
    return ( i == 0 ) ? I : ( ( i == 1 ) ? J : 1 );
  }

  /// Return value size
  static inline auto value_size() -> size_t
  {
    return SIZE;
  }

private:
  static constexpr size_t SIZE { I * J };
  static constexpr size_t RANK { ( I > 1 ) ? ( J > 1 ? 2 : 1 ) : 0 };
};

//-----------------------------------------------------------------------------

} // end namespace dolfin

#endif /* __DOLFIN_FUNCTION_VALUE_SPACE_H_ */
