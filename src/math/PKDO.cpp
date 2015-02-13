// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-06-03
// Last changed: 2008-04-22

#include <dolfin/math/PKDO.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/math/Jacobi.h>

#include <cmath>

namespace dolfin
{

//-----------------------------------------------------------------------------
PKDO::PKDO(uint i, uint j) :
    i_(i),
    j_(j)
{
}
//-----------------------------------------------------------------------------
real PKDO::eval(uint i, uint j, real r, real s)
{
  dolfin_assert(r >= 0.0);
  dolfin_assert(s >= -1.0);
  dolfin_assert(r + s <= 0.0);
  real cij = std::sqrt(0.5 * (2 * i + 1) * (i + j + 1));
  real rQ = (2 * r + s + 1) / (1 - s);
  return cij * Jacobi::eval(i, 0, 0, rQ) * std::pow(0.5 * (1 - s), (real) i)
      * Jacobi::eval(j, 2 * i + 1, 0, s);
}
//-----------------------------------------------------------------------------
real PKDO::ddx(uint i, uint j, real r, real s)
{
  error("PKDO derivatives not implemented.");
  return 0.0;
}
//-----------------------------------------------------------------------------

}
