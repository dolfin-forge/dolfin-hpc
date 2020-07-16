// Copyright (C) 2003 Johan Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_ARRAY_H
#define __DOLFIN_ARRAY_H

#include <dolfin/common/types.h>

#include <algorithm>
#include <vector>

namespace dolfin
{

template< typename T >
using Array = std::vector< T >;

//-----------------------------------------------------------------------------

template< typename T >
void destruct( Array< T * > & objects )
{
  for ( T * element : objects )
    if ( element != nullptr )
      delete element;

  objects.clear();
}

//-----------------------------------------------------------------------------

template < typename T, typename Iterator >
inline void append( Array< T > & array, Iterator begin, Iterator end )
{
#ifdef __SUNPRO_CC
 for ( Iterator it = begin; it != end; ++it )
 {
   array.push_back( *it );
 }
#else
 array.insert( array.end(), begin, end );
#endif
}

//-----------------------------------------------------------------------------

template < typename T >
inline uint max_array_size( Array< Array< T > > const & arrays )
{
    return std::max_element( arrays.begin(), arrays.end(),
                             []( Array< T > const & a, Array< T > const & b ) {
                               return a.size() < b.size();
                             } )->size();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_ARRAY_H */
