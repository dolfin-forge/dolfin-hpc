// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
// This file provides DOLFIN typedefs for basic types.

#ifndef __DOLFIN_TYPES_H
#define __DOLFIN_TYPES_H

#include <dolfin/config/dolfin_config.h>

#if HAVE_PARALLEL_HASH_MAP
#include <parallel_hashmap/phmap.h>
#elif __sgi
#include <hash_map>
#include <hash_set>
#else
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
using real    = double;

// Unsigned integers
using uint    = std::size_t;
using size_t  = std::size_t;

// Index type (at least 64bit)
using uidx    = uint64_t;

// Complex numbers
using complex = std::complex< double >;

//-----------------------------------------------------------------------------

constexpr uint DOLFIN_UINT_MIN   = std::numeric_limits< uint >::min();
constexpr uint DOLFIN_UINT_MAX   = std::numeric_limits< uint >::max();
constexpr uint DOLFIN_UINT_UNDEF = std::numeric_limits< uint >::max();

constexpr int DOLFIN_INT_MIN   = std::numeric_limits< int >::min();
constexpr int DOLFIN_INT_MAX   = std::numeric_limits< int >::max();
constexpr int DOLFIN_INT_UNDEF = std::numeric_limits< int >::max();

constexpr real DOLFIN_REAL_MIN   = std::numeric_limits< real >::min();
constexpr real DOLFIN_REAL_MAX   = std::numeric_limits< real >::max();
constexpr real DOLFIN_REAL_UNDEF = std::numeric_limits< real >::max();

constexpr long DOLFIN_LONG_MIN   = std::numeric_limits< long >::min();
constexpr long DOLFIN_LONG_MAX   = std::numeric_limits< long >::max();
constexpr long DOLFIN_LONG_UNDEF = std::numeric_limits< long >::max();

//-----------------------------------------------------------------------------

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

#if HAVE_PARALLEL_HASH_MAP

template < typename Key,
           typename Value,
           typename Hash  = phmap::priv::hash_default_hash< Key >,
           typename Eq    = phmap::priv::hash_default_eq< Key >,
           typename Alloc = std::allocator< std::pair< const Key, Value > > >
using _map = phmap::flat_hash_map< Key, Value, Hash, Eq, Alloc >;

template < typename Key,
           typename Hash  = phmap::priv::hash_default_hash< Key >,
           typename Eq    = phmap::priv::hash_default_eq< Key >,
           typename Alloc = std::allocator< Key > >
using _set = phmap::flat_hash_set< Key, Hash, Eq, Alloc >;

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

#else

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

#endif

//-----------------------------------------------------------------------------

/// Facility to compare object through pointers
template < class T >
auto objptrcmp( T const * p0, T const * p1 ) -> bool
{
  if ( p0 == p1 )
  {
    return true;
  }
  else if ( ( p0 == nullptr && p1 != nullptr ) || ( p0 != nullptr && p1 == nullptr ) )
  {
    return false;
  }
  return ( *p0 == *p1 );
}

} // end namespace dolfin

#endif
