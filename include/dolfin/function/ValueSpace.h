#ifndef __LICORNE_FUNCTION_VALUE_SPACE_H_
#define __LICORNE_FUNCTION_VALUE_SPACE_H_

#include <dolfin/common/types.h>

namespace dolfin
{

template<uint I = 1, uint J = 1>
class ValueSpace
{

  static uint const SIZE = I * J;

public:

  ValueSpace() :
    RANK((I > 1 ? (J > 1 ? 2 : 1) : 0))
  {
    DIM[0] = (RANK > 0 ? I : 1);
    DIM[1] = (RANK > 1 ? J : 1);
  }

  /// Return the rank of the value space
  uint rank() const
  {
    return RANK;
  }

  /// Return the dimension of the value space for axis i
  uint dim(uint i) const
  {
    return (i < RANK ? DIM[i] : 1);
  }

  /// Return value size
  uint value_size() const
  {
    return SIZE;
  }

private:

  uint DIM[2];
  uint const RANK;

};

} /* namespace licorne */

#endif /* __LICORNE_FUNCTION_VALUE_SPACE_H_ */
