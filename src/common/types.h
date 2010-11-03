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

namespace dolfin
{

  // Real numbers
  typedef double real;

  // Unsigned integers
  typedef unsigned int uint;

  // Complex numbers
  typedef std::complex<double> complex;

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
}

#endif
