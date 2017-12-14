// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-04-22
// Last changed: 2008-04-22
//
// This file provides DOLFIN typedefs for basic types.

#ifndef __DOLFIN_TYPES_H
#define __DOLFIN_TYPES_H

#include <dolfin/config/dolfin_config.h>

#include <complex>
#if (HAVE_TR1_UNORDERED_MAP && HAVE_TR1_UNORDERED_SET)
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#elif (__IBMCPP__ && __IBMCPP_TR1__)
#include <unordered_map>
#include <unordered_set>
#elif __sgi
#include <hash_map>
#include <hash_set>
#elif ENABLE_BOOST_TR1
#include <tr1/unordered_map.hpp>
#include <tr1/unordered_set.hpp>
#else
#include <map>
#include <set>
#endif

#include <climits>
#include <cfloat>
#include <stdint.h>

namespace dolfin
{

  // Real numbers
  typedef double real;

  // Unsigned integers
  typedef unsigned int uint;

  // Index type (at least 64bit)
  typedef uint64_t uidx;

  // Complex numbers
  typedef std::complex<double> complex;

  uint const DOLFIN_UINT_MIN = 0;
  uint const DOLFIN_UINT_MAX = UINT_MAX;
  uint const DOLFIN_UINT_UNDEF = UINT_MAX;

  int const DOLFIN_INT_MIN = INT_MIN;
  int const DOLFIN_INT_MAX = INT_MAX;
  int const DOLFIN_INT_UNDEF = INT_MAX;

  real const DOLFIN_REAL_MIN = DBL_MIN;
  real const DOLFIN_REAL_MAX = DBL_MAX;
  real const DOLFIN_REAL_UNDEF = DBL_MAX;

#if (HAVE_TR1_UNORDERED_MAP && HAVE_TR1_UNORDERED_SET)
#define _map std::tr1::unordered_map
#define _set std::tr1::unordered_set
#elif (__IBMCPP__ && __IBMCPP_TR1__)
#define _map std::tr1::unordered_map
#define _set std::tr1::unordered_set
#elif __sgi
#define _map std::hash_map
#define _set std::hash_set
#elif ENABLE_BOOST_TR1
#define _map std::tr1::unordered_map
#define _set std::tr1::unordered_set
#else
#define _map std::map
#define _set std::set
#endif

//-----------------------------------------------------------------------------

/// Facility to compare arrays
template <class T> bool cmp(size_t N, T const * x0, T const * x1)
{
  if (x0 == x1)
  {
    return true;
  }
  else if ((x0 == NULL || x1 == NULL))
  {
    return false;
  }
  for (size_t ii = 0; ii < N; ++ii)
  {
    if (x0[ii] != x1[ii])
    {
      return false;
    }
  }
  return true;
}

/// Facility to compare arrays
template <class T> bool cmp(T const * x0b, T const * x0e, T const * x1)
{
  return cmp(x0e - x0b, x0b, x1);
}

/// Facility to compare object through pointers
template<class T> bool objptrcmp(T const * p0, T const * p1)
{
  if (p0 == p1)
  {
    return true;
  }
  else if ((p0 == NULL && p1 != NULL) || (p0 != NULL && p1 == NULL))
  {
    return false;
  }
  return (*p0 == *p1);
}

/// Equivalence operator for unordered maps to be used for assertion checking
template<class K, class V> bool operator==(_map<K,V> const& m0,
                                           _map<K,V> const& m1)
{
  // The distance from begin to end should be the same
  if (m0.size() != m1.size())
  {
    return false;
  }
  // Both groups returned by equal range have equal size and there exists
  // a permutation such that the elements are equal two by two.
  // In this case this is just pair comparison.
  typename _map<K, V>::const_iterator it1;
  for (typename _map<K, V>::const_iterator it0 = m0.begin(); it0 != m0.end();
       ++it0)
  {
    it1 = m1.find(it0->first);
    if(it1 != m1.end())
    {
      if(it1->second != it0->second)
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }
  // Alrighty
  return true;
}

/// !Equivalence operator for unordered maps to be used for assertion checking
template<class K, class V> bool operator!=(_map<K,V> const& m0,
                                           _map<K,V> const& m1)
{
  return !(m0 == m1);
}

/// Equivalence operator for unordered maps to be used for assertion checking
template<class T> bool operator==(_set<T> const& m0, _set<T> const& m1)
{
  // The distance from begin to end should be the same
  if (m0.size() != m1.size())
  {
    return false;
  }
  // Both groups returned by equal range have equal size and there exists
  // a permutation such that the elements are equal two by two.
  // In this case this is just an element comparison.
  typename _set<T>::const_iterator it1;
  for (typename _set<T>::const_iterator it0 = m0.begin(); it0 != m0.end();
       ++it0)
  {
    it1 = m1.find(*it0);
    if(it1 != m1.end())
    {
      if(*it1 != *it0)
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }
  // Alrighty
  return true;
}

/// !Equivalence operator for unordered maps to be used for assertion checking
template<class T> bool operator!=(_set<T> const& m0, _set<T> const& m1)
{
  return !(m0 == m1);
}

/// Unordered set intersection
template<class T>
void intersection(_set<T> const& s0, _set<T> const& s1, _set<T>& in)
{
  in = s0;
  for(_set<uint>::iterator it = in.begin(); it != in.end();)
  {
    if (s1.count(*it) == 0)
    {
      in.erase(it++);
    }
    else
    {
      ++it;
    }
  }
}

}

#endif
