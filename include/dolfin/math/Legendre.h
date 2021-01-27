// Copyright (C) 2003-2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LEGENDRE_H
#define __DOLFIN_LEGENDRE_H

#include <dolfin/common/types.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

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
  Legendre( size_t n );

  ///
  ~Legendre() = default;

  /// Evaluation at given point
  auto operator()( real x ) -> real;

  /// Evaluation of derivative at given point
  auto ddx( real x ) -> real;

  /// Evaluation of second derivative at given point
  auto d2dx( real x ) -> real;

  /// Evaluation at given point
  static auto eval( size_t n, real x ) -> real;

  /// Evaluation of derivative at given point
  static auto ddx( size_t n, real x ) -> real;

  /// Evaluation of second derivative at given point
  static auto d2dx( size_t n, real x ) -> real;

private:
  size_t const n_;
};

//--- INLINES -----------------------------------------------------------------

inline auto Legendre::operator()( real x ) -> real
{
  return eval( n_, x );
}

//-----------------------------------------------------------------------------

inline auto Legendre::ddx( real x ) -> real
{
  return ddx( n_, x );
}

//-----------------------------------------------------------------------------

inline auto Legendre::d2dx( real x ) -> real
{
  return d2dx( n_, x );
}

//-----------------------------------------------------------------------------

} // namespace dolfin

#endif
