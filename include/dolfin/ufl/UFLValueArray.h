// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_VALUE_ARRAY_H
#define __UFL_VALUE_ARRAY_H

#include <string>
#include <vector>

#include <dolfin/common/types.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  ValueArray
 *
 *  @brief  Provides a basic type for arrays of natural numbers used for value
 *          shapes.
 */

class ValueArray : public std::vector<dolfin::uint>
{

public:

  ///
  ValueArray();

  ///
  ValueArray(dolfin::uint const i);

  ///
  ValueArray(dolfin::uint const k, dolfin::uint const i);

  ///
  ~ValueArray();

  /// Return the product of all entries
  dolfin::uint prod() const;

  ///
  std::string str() const;

private:

};

/// Returns copy appending second value array to the first
ValueArray operator+(ValueArray const&v1, ValueArray const& v2);

} /* namespace ufl */

#endif /* __UFL_VALUE_ARRAY_H */

