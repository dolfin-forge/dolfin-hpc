// Copyright (C) 2003-2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_BASIC_H
#define __DOLFIN_BASIC_H

#include <dolfin/common/constants.h>
#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>

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
static inline bool abscmp(real x, real y)
{
  return std::fabs(x - y) < DOLFIN_EPS;
}

/// Convenient function
static inline bool small(real x, real eps = DOLFIN_EPS)
{
  return std::fabs(x) < eps;
}

/// Return absolute real comparison for ~ O(1) with given epsilon:
/// | x - y | < eps
static inline bool abscmp(real x, real y, real eps)
{
  return std::fabs(x - y) < std::fabs(eps);
}

/// Return weak relative real comparison:
/// ( |(x - y) / x| < eps ) || ( |(x - y) / y| < eps )
static inline bool wrelcmp(real x, real y, real eps)
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
static inline bool srelcmp(real x, real y, real eps)
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
static inline real sqr(real x)
{
  return x * x;
}

/// Return the euclidian norm
static inline real norm1(uint n, real * x)
{
  real val = 0.0;
  for (uint i = 0; i < n; ++i)
  {
    val += std::fabs(x[i]);
  }
  return val;
}

/// Return the euclidian norm
static inline real norm2(uint n, real * x)
{
  real val = 0.0;
  for (uint i = 0; i < n; ++i)
  {
    val += x[i]*x[i];
  }
  return std::sqrt(val);
}

/// Return a to the power n
static inline uint ipow(uint a, uint n)
{
  uint p = a;
  for (uint i = 1; i < n; i++)
    p *= a;
  return p;
}

/// Return factorial of a
static inline int fact(int a, int acc = 1)
{
  // error condition
  if (acc < 0)
  {
    return -1;
  }

  // termination condition
  if (a == 0||a == 1)
    return acc;

  // Tail recursive call
  return fact(a - 1, acc * a);
}

/// Gamma function from PELICANS
static inline real gamma(double const x)
{
  int n = x < 1.5 ? -((int) (2.5 - x)) : (int) (x - 1.5);
  double w = x - (n + 2);
  double y = ((((((((((((-1.99542863674e-7 * w + 1.337767384067e-6) * w
      - 2.591225267689e-6) * w - 1.7545539395205e-5) * w + 1.45596568617526e-4)
      * w - 3.60837876648255e-4) * w - 8.04329819255744e-4) * w
      + 0.008023273027855346) * w - 0.017645244547851414) * w
      - 0.024552490005641278) * w + 0.19109110138763841) * w
      - 0.233093736421782878) * w - 0.422784335098466784) * w
      + 0.99999999999999999;
  if (n > 0)
  {
    w = x - 1;
    for (int k = 2; k <= n; k++)
    {
      w *= x - k;
    }
  }
  else
  {
    w = 1;
    for (int k = 0; k > n; k--)
    {
      y *= x - k;
    }
  }
  double result = w / y;
  return (result);
}

/// Seed only first time
static bool rand_seeded = false;

/// Return a random number, uniformly distributed between [0.0, 1.0)
/// !!! Not quite, the implementation does not ensure that strongly.
static inline real rand()
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
static inline void seed(unsigned int s)
{
  std::srand(s);
  rand_seeded = true;
}

///
static inline real percent(uint n, uint d)
{
  return (100.0 * real(n) / real(d));
}

/// Return sequence of natural numbers ranging from begin to (excluding) end
template<class Iterator, class T>
static inline void range(Iterator begin, Iterator end, T v = T(), int s = 1)
{
  if (s)
  {
    if (s > 0)
    {
      while (begin != end)
      {
        *begin++ = v;
        for (int i = 0; i < s; ++i) ++v;
      }
    }
    else
    {
      while (begin != end)
      {
        *begin++ = v;
        for (int i = 0; i > s; --i) --v;
      }
    }
  }
  else
  {
    std::fill(begin, end, v);
  }
}

/// e0 contains e1
static inline bool contains( Array<uint> const & e0, Array<uint> const & e1 )
{
  for ( Array<uint>::const_iterator e = e1.begin(); e != e1.end(); ++e )
    if ( e0.end() == std::find( e0.begin(), e0.end(), *e ) )
      return false;
  return true;
}

/// Return rank within linear distribution
static inline uint rank(uint L, uint R, uint i)
{
  return static_cast<uint>(std::max(
      std::floor((double) i / (double) (L + 1)),
      std::floor((double) ((double) i - (double) R) / (double) L)));
}

/// Return the number of significant digits in the integer
static inline uint sdigits(uint n)
{
  uint i = 0;
  while (n > 0)
  {
    ++i;
    n /= 10;
  }
  return i;
}

/// Check if Not-a-Number if C99 is not used
#if __STDC_VERSION__ < 199901L
#ifndef isnan
static inline bool isnan(real x)
{
  return x != x;
}
#endif
#endif

}

#endif
