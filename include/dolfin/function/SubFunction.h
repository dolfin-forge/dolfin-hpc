// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-04-27
// Last changed: 2007-04-29

#ifndef __SUB_FUNCTION_H
#define __SUB_FUNCTION_H

#include <dolfin/common/types.h>

namespace dolfin
{

class DiscreteFunction;

/// This class represents a sub function (view) of a (discrete function).
/// It's purpose is to enable expressions like
///
///    Function w;
///    Function u = w[0];
///    Function p = w[1];
///
/// without needing to create and destroy temporaries. No data is created
/// until a Function is assigned to a SubFunction, at which point the data
/// needed to represent the sub function is created.

class SubFunction
{
public:

  /// Create sub function
  SubFunction(DiscreteFunction& f, uint i) :
      f_(&f),
      i_(i)
  {
  }

  /// Destructor
  ~SubFunction()
  {
  }

  /// Return global function
  DiscreteFunction& function() const;

  /// Return index of the sub function
  uint index() const;

  /// Display basic information
  void disp() const;

private:

  /// Create empty sub function
  SubFunction() :
      f_(NULL),
      i_(0)
  {
  }

  // Pointer to discrete function
  DiscreteFunction * const f_;

  // Sub function index
  uint const i_;

};

//-----------------------------------------------------------------------------
inline DiscreteFunction& SubFunction::function() const
{
  return *f_;
}

//-----------------------------------------------------------------------------
inline uint SubFunction::index() const
{
  return i_;
}

}

#endif
