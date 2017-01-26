// Copyright (C) 2013 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-05-24 (merged from branch larcher)
// Last changed: 2013-05-24

#ifndef __DOLFIN_EXPRESSION_H
#define __DOLFIN_EXPRESSION_H

#include <dolfin/common/types.h>

namespace dolfin
{

/**
 *  @class  Expression
 *
 *  @brief  Derived classes should implement the copy constructor if they
 *          allocate new objects.
 */

class Expression
{
public:

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const = 0;

  /// Return the rank of the value space
  virtual uint rank() const = 0;

  /// Return the dimension of the value space for axis i
  virtual uint dim(uint i) const = 0;

  // Return the value size
  inline uint value_size() const
  {
    uint size = 1;
    for (uint i = 0; i < this->rank(); ++i)
    {
      size *= this->dim(i);
    }
    return size;
  }

protected:

  /// Create user-defined function
  Expression() {}

  /// Destructor
  virtual ~Expression() {}

};

} /* namespace licorne */

#endif /* __DOLFIN_EXPRESSION_H */
