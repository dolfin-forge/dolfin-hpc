// Copyright (C) 2003-2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-02-06
// Last changed: 2006-06-16

#ifndef __BASIC_H
#define __BASIC_H

#include <dolfin/common/types.h>

#include <time.h>
#include <cstdlib>
#include <cmath>

namespace dolfin
{

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
/// !!! Not quite
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
