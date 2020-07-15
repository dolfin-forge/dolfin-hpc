// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
// This file provides DOLFIN typedefs for basic types.

#ifndef __DOLFIN_TYPES_H
#define __DOLFIN_TYPES_H

#include <dolfin/config/dolfin_config.h>

#if HAVE_PARALLEL_HASH_MAP
#include <parallel_hashmap/phmap.h>
#elif ( HAVE_TR1_UNORDERED_MAP && HAVE_TR1_UNORDERED_SET )
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#elif ( __IBMCPP__ && __IBMCPP_TR1__ )
#include <unordered_map>
#include <unordered_set>
#elif __sgi
#include <hash_map>
#include <hash_set>
#elif ( HAVE_UNORDERED_MAP && HAVE_UNORDERED_SET )
#include <unordered_map>
#include <unordered_set>
#endif

#include <cfloat>
#include <complex>
#include <cstdint>
#include <limits>
#include <map>
#include <set>

namespace dolfin
{

// Real numbers
typedef double real;

// Unsigned integers
typedef unsigned int uint;

// Index type (at least 64bit)
typedef uint64_t uidx;

// Complex numbers
typedef std::complex< double > complex;

//-----------------------------------------------------------------------------

uint const DOLFIN_UINT_MIN   = std::numeric_limits< uint >::min();
uint const DOLFIN_UINT_MAX   = std::numeric_limits< uint >::max();
uint const DOLFIN_UINT_UNDEF = std::numeric_limits< uint >::max();

int const DOLFIN_INT_MIN   = std::numeric_limits< int >::min();
int const DOLFIN_INT_MAX   = std::numeric_limits< int >::max();
int const DOLFIN_INT_UNDEF = std::numeric_limits< int >::max();

real const DOLFIN_REAL_MIN   = std::numeric_limits< real >::min();
real const DOLFIN_REAL_MAX   = std::numeric_limits< real >::max();
real const DOLFIN_REAL_UNDEF = std::numeric_limits< real >::max();

long const DOLFIN_LONG_MIN   = std::numeric_limits< long >::min();
long const DOLFIN_LONG_MAX   = std::numeric_limits< long >::max();
long const DOLFIN_LONG_UNDEF = std::numeric_limits< long >::max();

//-----------------------------------------------------------------------------

#if HAVE_PARALLEL_HASH_MAP

template < typename Key,
           typename Value,
           typename Hash  = phmap::container_internal::hash_default_hash< Key >,
           typename Eq    = phmap::container_internal::hash_default_eq< Key >,
           typename Alloc = std::allocator< std::pair< const Key, Value > > >
using _map = phmap::flat_hash_map< Key, Value, Hash, Eq, Alloc >;

template < typename Key,
           typename Hash  = phmap::container_internal::hash_default_hash< Key >,
           typename Eq    = phmap::container_internal::hash_default_eq< Key >,
           typename Alloc = std::allocator< Key > >
using _set = phmap::flat_hash_set< Key, Hash, Eq, Alloc >;

template < typename Key,
           typename Value,
           typename Compare = std::less< Key >,
           typename Allocator =
             std::allocator< std::pair< const Key, Value > > >
using _ordered_map = std::map< Key, Value, Compare, Allocator >;

template < typename Key,
           typename Compare   = std::less< Key >,
           typename Allocator = std::allocator< Key > >
using _ordered_set = std::set< Key, Compare, Allocator >;

#else

template < typename Key,
           typename Value,
           typename Comp  = std::less< Key >,
           typename Alloc = std::allocator< std::pair< const Key, Value > > >
using _ordered_map = std::map< Key, Value, Comp, Alloc >;

template < typename Key,
           typename Comp  = std::less< Key >,
           typename Alloc = std::allocator< Key > >
using _ordered_set = std::set< Key, Comp, Alloc >;

#if ( HAVE_TR1_UNORDERED_MAP && HAVE_TR1_UNORDERED_SET ) \
  || ( __IBMCPP__ && __IBMCPP_TR1__ )

template < typename Key,
           typename Value,
           typename Hash  = std::hash< Key >,
           typename Comp  = std::equal_to< Key >,
           typename Alloc = std::allocator< std::pair< const Key, Value > > >
using _map = std::tr1::unordered_map< Key, Value, Hash, Comp, Alloc >;

template < typename Key,
           typename Hash  = std::hash< Key >,
           typename Comp  = std::equal_to< Key >,
           typename Alloc = std::allocator< Key > >
using _set = std::tr1::unordered_set< Key, Hash, Comp, Alloc >;

#elif __sgi

template < typename Key,
           typename Value,
           typename Hash  = std::hash< Key >,
           typename Comp  = std::equal_to< Key >,
           typename Alloc = std::allocator< std::pair< const Key, Value > > >
using _map = std::hash_map< Key, Value, Hash, Comp, Alloc >;

template < typename Key,
           typename Hash  = std::hash< Key >,
           typename Comp  = std::equal_to< Key >,
           typename Alloc = std::allocator< Key > >
using _set = std::hash_set< Key, Hash, Comp, Alloc >;

#elif ( HAVE_UNORDERED_MAP && HAVE_UNORDERED_SET )

template < typename Key,
           typename Value,
           typename Hash  = std::hash< Key >,
           typename Comp  = std::equal_to< Key >,
           typename Alloc = std::allocator< std::pair< const Key, Value > > >
using _map = std::unordered_map< Key, Value, Hash, Comp, Alloc >;

template < typename Key,
           typename Hash  = std::hash< Key >,
           typename Comp  = std::equal_to< Key >,
           typename Alloc = std::allocator< Key > >
using _set = std::unordered_set< Key, Hash, Comp, Alloc >;

#else

template < typename Key,
           typename Value,
           typename Comp  = std::less< Key >,
           typename Alloc = std::allocator< std::pair< const Key, Value > > >
using _map = std::map< Key, Value, Comp, Alloc >;

template < typename Key,
           typename Comp  = std::less< Key >,
           typename Alloc = std::allocator< Key > >
using _set = std::set< Key, Comp, Alloc >;

#endif

#endif

//-----------------------------------------------------------------------------

/// Facility to compare object through pointers
template < class T >
bool objptrcmp( T const * p0, T const * p1 )
{
  if ( p0 == p1 )
  {
    return true;
  }
  else if ( ( p0 == NULL && p1 != NULL ) || ( p0 != NULL && p1 == NULL ) )
  {
    return false;
  }
  return ( *p0 == *p1 );
}

} // end namespace dolfin

#endif
