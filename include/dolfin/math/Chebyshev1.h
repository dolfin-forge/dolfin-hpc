// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CHEBYSHEV_FIRST_H
#define __DOLFIN_CHEBYSHEV_FIRST_H

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/types.h>
#include <dolfin/math/Jacobi.h>

namespace dolfin
{

/// Chebyshev polynomials of the first kind

class Chebyshev1
{

public:

  /// Constructor
  Chebyshev1(uint n) :
      n_(n)
  {
  }

  /// Destructor
  ~Chebyshev1()
   {
   }

  /// Evaluation at given point
  real operator()(real x);

  /// Evaluation of derivative at given point
  real ddx(real x);

  /// Evaluation at given point
  static real eval(uint n, real x);

  /// Evaluation of derivative at given point
  static real ddx(uint n, real x);

private:

  uint const n_;

};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline real Chebyshev1::operator()(real x)
{
  return eval(n_, x);
}
//-----------------------------------------------------------------------------
inline real Chebyshev1::ddx(real x)
{
  return ddx(n_, x);
}

//-----------------------------------------------------------------------------
inline real eval(uint n, real x)
{
  return Jacobi::eval(n, -0.5, -0.5, x);
}

//-----------------------------------------------------------------------------
inline real ddx(uint n, real x)
{
  return Jacobi::ddx(n, -0.5, -0.5, x);
}

}

#endif
