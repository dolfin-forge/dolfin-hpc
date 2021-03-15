// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CHEBYSHEV_SECOND_H
#define __DOLFIN_CHEBYSHEV_SECOND_H

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/types.h>
#include <dolfin/math/Jacobi.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

/// Chebyshev polynomials of the second kind

class Chebyshev2
{

public:

  /// Constructor
  Chebyshev2(size_t n) :
      n_(n)
  {
  }

  /// Destructor
  ~Chebyshev2()
  {
  }

  /// Evaluation at given point
  real operator()(real x);

  /// Evaluation of derivative at given point
  real ddx(real x);

  /// Evaluation at given point
  static real eval(size_t n, real x);

  /// Evaluation of derivative at given point
  static real ddx(size_t n, real x);

private:

  size_t const n_;

};

//--- INLINES -----------------------------------------------------------------

inline real Chebyshev2::operator()(real x)
{
  return eval(n_, x);
}

//-----------------------------------------------------------------------------

inline real Chebyshev2::ddx(real x)
{
  return ddx(n_, x);
}

//-----------------------------------------------------------------------------

inline real eval(size_t n, real x)
{
  return Jacobi::eval(n, +0.5, +0.5, x);
}

//-----------------------------------------------------------------------------

inline real ddx(size_t n, real x)
{
  return Jacobi::ddx(n, +0.5, +0.5, x);
}

//-----------------------------------------------------------------------------

} // namespace dolfin

#endif
