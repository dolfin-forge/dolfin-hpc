// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/math/Jacobi.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>

#include <cmath>

namespace dolfin
{

//-----------------------------------------------------------------------------
Jacobi::Jacobi(uint n, real a, real b) :
    n_(n),
    a_(a),
    b_(b)
{
}
//-----------------------------------------------------------------------------
real Jacobi::eval(uint n, real a, real b, real x)
{
  // Special case n = 0
  if (n == 0) return 1.0;

  // Special case n = 1
  if (n == 1) return 0.5 * (a - b + (a + b + 2.0) * x);

  //
  real p = 0.0;
  real const apb = a + b;
  real pn2 = 1.0;
  real pn1 = 0.5 * (a - b + (apb + 2.0) * x);
  for (uint k = 0; k < n + 1; ++k)
  {
    real a1 = 2.0 * k * (k + apb) * (2.0 * k + apb - 2.0);
    real a2 = (2.0 * k + apb - 1.0) * (a * a - b * b);
    real a3 = (2.0 * k + apb - 2.0) * (2.0 * k + apb - 1.0) * (2.0 * k + apb);
    real a4 = 2.0 * (k + a - 1.0) * (k + b - 1.0) * (2.0 * k + apb);
    a2 = a2 / a1;
    a3 = a3 / a1;
    a4 = a4 / a1;
    p = (a2 + a3 * x) * pn1 - a4 * pn2;
    pn2 = pn1;
    pn1 = p;
  }
  return p;
}
//-----------------------------------------------------------------------------
real Jacobi::ddx(uint n, real a, real b, real x)
{
  // Special case n = 0
  if (n == 0) return 0.0;

  //
  real nn = real(n);
  return 0.5 * (a + b + nn + 1) * eval(a + 1, b + 1, nn - 1, x);
}
//-----------------------------------------------------------------------------

}
