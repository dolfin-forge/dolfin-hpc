// Copyright (C) 2013 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.

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
  virtual void eval( real * values, const real * x ) const = 0;

  /// Return the rank of the value space
  virtual auto rank() const -> size_t = 0;

  /// Return the dimension of the value space for axis i
  virtual auto dim( size_t i ) const -> size_t = 0;

  // Return the value size
  virtual auto value_size() const -> size_t = 0;

protected:
  /// Create user-defined function
  Expression() = default;

  /// Destructor
  virtual ~Expression() = default;
};

} // end namespace dolfin

#endif /* __DOLFIN_EXPRESSION_H */
