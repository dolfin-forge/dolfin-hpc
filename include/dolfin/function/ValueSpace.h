#ifndef __LICORNE_FUNCTION_VALUE_SPACE_H_
#define __LICORNE_FUNCTION_VALUE_SPACE_H_

#include <dolfin/common/types.h>

namespace dolfin
{

template<uint I = 1, uint J = 1>
class ValueSpace
{

public:

  ValueSpace()
  {
  }

  /// Return the rank of the value space
  static inline uint rank()
  {
    return RANK;
  }

  /// Return the dimension of the value space for axis i
  static inline uint dim(uint i)
  {
    return (i < RANK ? DIM[i] : 1);
  }

  /// Return value size
  static inline uint value_size()
  {
    return SIZE;
  }

private:

  static uint const SIZE = I * J;
  static uint const DIM[2];
  static uint const RANK;

};

template<uint I, uint J>
uint const ValueSpace<I, J>::DIM[2] = { I, J };

template<uint I, uint J>
uint const ValueSpace<I, J>::RANK   = (I > 1 ? (J > 1 ? 2 : 1) : 0);

} // end namespace dolfin

#endif /* __LICORNE_FUNCTION_VALUE_SPACE_H_ */
