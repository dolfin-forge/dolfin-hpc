// Copyright (C) 2003-2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-02-06
// Last changed: 2006-06-16

#ifndef __BASIC_H
#define __BASIC_H

#include <dolfin/common/constants.h>
#include <dolfin/common/types.h>

#include <time.h>
#include <cstdlib>
#include <cmath>
#include <limits>

namespace dolfin
{

/// References:
/// [1] Knuth D.E. The art of computer programming (vol II).
/// [2] Boost: The Test Tools, Copyright Gennadiy Rozental 2001-2005.

/// Return absolute real comparison for ~ O(1) with default epsilon:
/// | x - y | < eps
inline bool abscmp(real x, real y)
{
  return std::fabs(x - y) < DOLFIN_EPS;
}

/// Return absolute real comparison for ~ O(1) with given epsilon:
/// | x - y | < eps
inline bool abscmp(real x, real y, real eps)
{
  return std::fabs(x - y) < std::fabs(eps);
}

/// Return weak relative real comparison:
/// ( |(x - y) / x| < eps ) || ( |(x - y) / y| < eps )
inline bool wrelcmp(real x, real y, real eps)
{
  real const d = std::fabs(x - y);
  real const m = std::max(std::fabs(x), std::fabs(y));
  // Take very small positive values of d and m into account
  if (d < std::numeric_limits<real>::min() ||
      m < std::numeric_limits<real>::min())
  {
    return true;
  }
  // Trying to avoid underflow issues in most common cases
  return (
      (m > DOLFIN_EPS) ?
          (d / m < std::fabs(eps)) :
          (d / std::sqrt(m) < std::fabs(eps) * std::sqrt(m)));
}

/// Return strong relative real comparison:
/// ( |(x - y) / x| < eps ) && ( |(x - y) / y| < eps )
inline bool srelcmp(real x, real y, real eps)
{
  real const d = std::fabs(x - y);
  real const m = std::min(std::fabs(x), std::fabs(y));
  // Take very small positive values of d and m into account
  if (d < std::numeric_limits<real>::min() ||
      m < std::numeric_limits<real>::min())
  {
    return true;
  }
  // Trying to avoid underflow issues in most common cases
  return (
      (m > DOLFIN_EPS) ?
          (d / m < std::fabs(eps)) :
          (d / std::sqrt(m) < std::fabs(eps) * std::sqrt(m)));
}

/// Return the square of x
inline real sqr(real x)
{
  return x * x;
}

/// Return a to the power n
inline uint ipow(uint a, uint n)
{
  uint p = a;
  for (uint i = 1; i < n; i++)
    p *= a;
  return p;
}

/// Seed only first time
static bool rand_seeded = false;

/// Return a random number, uniformly distributed between [0.0, 1.0)
/// !!! Not quite, the implementation does not ensure that strongly.
inline real rand()
{
  if (!rand_seeded)
  {
    unsigned int s = static_cast<long int>(::time(0));
    std::srand(s);
    rand_seeded = true;
  }

  return static_cast<real>(std::rand()) / static_cast<real>(RAND_MAX);
}

/// Seed random number generator
inline void seed(unsigned int s)
{
  std::srand(s);
  rand_seeded = true;
}

}

#endif
