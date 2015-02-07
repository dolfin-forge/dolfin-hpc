// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-11-08
// Last changed: 2014-12-08

#ifndef __JACOBI_H
#define __JACOBI_H

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/types.h>

namespace dolfin
{
/// Jacobi polynomial from FIAT 1.4.0

class Jacobi
{

public:

  Jacobi(uint n, real a, real b);

  /// Evaluation at given point
  real operator()(real x);

  /// Evaluation of derivative at given point
  real ddx(real x);

  /// Evaluation at given point with weights a , b
  static real eval(uint n, real a, real b, real x);

  /// Evaluation of derivative at given point with weights a , b
  static real ddx(uint n, real a, real b, real x);

private:

  uint const n_;
  real const a_;
  real const b_;

};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline real Jacobi::operator() (real x)
{
  return eval(n_, a_, b_, x);
}
//-----------------------------------------------------------------------------
inline real Jacobi::ddx(real x)
{
  return ddx(n_, a_, b_, x);
}

}

#endif
