// Copyright (C) 2003-2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LEGENDRE_H
#define __DOLFIN_LEGENDRE_H

#include <dolfin/common/types.h>

namespace dolfin
{

/// Legendre polynomial of given degree n on the interval [-1,1].
///
///   P0(x) = 1
///   P1(x) = x
///   P2(x) = (3x^2 - 1) / 2
///   ...
///
/// The function values and derivatives are computed using
/// three-term recurrence formulas.

class Legendre
{

public:

  ///
  Legendre(uint n);

  ///
  ~Legendre();

  /// Evaluation at given point
  real operator()(real x);

  /// Evaluation of derivative at given point
  real ddx(real x);

  /// Evaluation of second derivative at given point
  real d2dx(real x);

  /// Evaluation at given point
  static real eval(uint n, real x);

  /// Evaluation of derivative at given point
  static real ddx(uint n, real x);

  /// Evaluation of second derivative at given point
  static real d2dx(uint n, real x);

private:

  uint const n_;

};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline real Legendre::operator()(real x)
{
  return eval(n_, x);
}
//-----------------------------------------------------------------------------
inline real Legendre::ddx(real x)
{
  return ddx(n_, x);
}
//-----------------------------------------------------------------------------
inline real Legendre::d2dx(real x)
{
  return d2dx(n_, x);
}

}

#endif
