// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/math/Legendre.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/math/Jacobi.h>

#include <cmath>

namespace dolfin
{

//-----------------------------------------------------------------------------
Legendre::Legendre(uint n) :
    n_(n)
{
}
//-----------------------------------------------------------------------------
auto Legendre::eval(uint n, real x) -> real
{
  // Special case n = 0
  if (n == 0) return 1.0;

  // Special case n = 1
  if (n == 1) return x;

  // Use Jacobi with a = b =0 to evaluate Legendre
  return Jacobi::eval(n, 0.0, 0.0, x);
}
//-----------------------------------------------------------------------------
auto Legendre::ddx(uint n, real x) -> real
{
  // Special case n = 0
  if (n == 0) return 0.0;

  // Special case n = 1
  if (n == 1) return 1.0;

  // Avoid division by zero
  if (fabs(x - 1.0) < DOLFIN_EPS) x -= 2.0 * DOLFIN_EPS;
  if (fabs(x + 1.0) < DOLFIN_EPS) x += 2.0 * DOLFIN_EPS;

  // Formula, BETA page 254
  real nn = real(n);
  return nn * (x * eval(n, x) - eval(n - 1, x)) / (x * x - 1.0);
}
//-----------------------------------------------------------------------------
auto Legendre::d2dx(uint n, real x) -> real
{
  // Special case n = 0
  if (n == 0) return 0.0;

  // Special case n = 1
  if (n == 1) return 0.0;

  // Avoid division by zero
  if (fabs(x - 1.0) < DOLFIN_EPS) x -= 2.0 * DOLFIN_EPS;
  if (fabs(x + 1.0) < DOLFIN_EPS) x += 2.0 * DOLFIN_EPS;

  // Formula, BETA page 254
  real nn = real(n);
  return (2.0 * x * ddx(n, x) - nn * (nn + 1) * eval(n, x)) / (1.0 - x * x);
}
//-----------------------------------------------------------------------------

}
