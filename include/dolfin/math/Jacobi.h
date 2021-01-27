// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_JACOBI_H
#define __DOLFIN_JACOBI_H

#include <dolfin/common/types.h>
#include <dolfin/log/dolfin_log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

/// Jacobi polynomial from FIAT 1.4.0

class Jacobi
{

public:
  Jacobi( size_t n, real a, real b );

  /// Evaluation at given point
  auto operator()( real x ) -> real;

  /// Evaluation of derivative at given point
  auto ddx( real x ) -> real;

  /// Evaluation at given point with weights a , b
  static auto eval( size_t n, real a, real b, real x ) -> real;

  /// Evaluation of derivative at given point with weights a , b
  static auto ddx( size_t n, real a, real b, real x ) -> real;

private:
  size_t const n_;
  real const   a_;
  real const   b_;
};

//--- INLINES -----------------------------------------------------------------

inline auto Jacobi::operator()( real x ) -> real
{
  return eval( n_, a_, b_, x );
}

//-----------------------------------------------------------------------------

inline auto Jacobi::ddx( real x ) -> real
{
  return ddx( n_, a_, b_, x );
}

//-----------------------------------------------------------------------------

} // namespace dolfin

#endif
