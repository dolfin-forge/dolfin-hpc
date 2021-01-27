#ifndef __LICORNE_FUNCTION_VALUE_SPACE_H_
#define __LICORNE_FUNCTION_VALUE_SPACE_H_

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
  static inline size_t rank()
  {
    return RANK;
  }

  /// Return the dimension of the value space for axis i
  static inline size_t dim( size_t i )
  {
    return ( i < RANK ? DIM[i] : 1 );
  }

  /// Return value size
  static inline size_t value_size()
  {
    return SIZE;
  }

private:
  static constexpr size_t SIZE { I * J };
  static constexpr size_t DIM[2] { I, J };
  static constexpr size_t RANK { ( I > 1 ) ? ( J > 1 ? 2 : 1 ) : 0 };
};

//-----------------------------------------------------------------------------

} // end namespace dolfin

#endif /* __LICORNE_FUNCTION_VALUE_SPACE_H_ */
