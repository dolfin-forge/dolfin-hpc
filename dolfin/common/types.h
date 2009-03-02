// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-04-22
// Last changed: 2008-04-22
//
// This file provides DOLFIN typedefs for basic types.

#ifndef __DOLFIN_TYPES_H
#define __DOLFIN_TYPES_H

#include <complex>
#ifdef __GNUG__
#if (__GNUG__ > 3 && __GXX_EXPERIMENTAL_CXX0X__) 
#include <unordered_map>
#include <unordered_set>
#else
#include <ext/hash_map>
#include <ext/hash_set>
#endif
#elif (__IBMCPP__ && __IBMCPP_TR1__) 
#include <unordered_map>
#include <unordered_set>
#elif __sgi
#include <hash_map>
#include <hash_set>
#else
#include <map>
#include <set>
#endif

namespace dolfin
{

  // Real numbers
  typedef double real;

  // Unsigned integers
  typedef unsigned int uint;

  // Complex numbers
  typedef std::complex<double> complex;

#ifdef __GNUG__
#if (__GNUG__ > 3 && __GXX_EXPERIMENTAL_CXX0X__)
#define _map std::unordered_map
#define _set std::unordered_set
#else
#define _map __gnu_cxx::hash_map 
#define _set __gnu_cxx::hash_set
#endif
#elif (__IBMCPP__ && __IBMCPP_TR1__) 
#define _map std::tr1::unordered_map
#define _set std::tr1::unordered_set
#elif __sgi
#define _map std::hash_map
#define _set std::hash_set
#else
#define _map std::map 
#define _set std::set
#endif
}

#endif
