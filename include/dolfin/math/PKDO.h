// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PKDO_H
#define __DOLFIN_PKDO_H

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/types.h>

namespace dolfin
{

/// Proriol-Koornwinder-Dubiner-Owens polynomials defined from the reference
/// triangle T = { (r,s), 0 >= r, s >= -1, r+ s <= 0 }

class PKDO
{

public:

  PKDO(uint i, uint j);

  /// Evaluation at given point
  auto operator()(real r, real s) -> real;

  /// Evaluation of derivative at given point
  auto ddx(real r, real s) -> real;

  /// Evaluation at given point with weights a , b
  static auto eval(uint i, uint j, real r, real s) -> real;

  /// Evaluation of derivative at given point with weights a , b
  static auto ddx(uint i, uint j, real r, real s) -> real;

private:

  uint const i_;
  uint const j_;

};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline auto PKDO::operator()(real r, real s) -> real
{
  return eval(i_, j_, r, s);
}
//-----------------------------------------------------------------------------
inline auto PKDO::ddx(real r, real s) -> real
{
  return ddx(i_, j_, r, s);
}

}

#endif
