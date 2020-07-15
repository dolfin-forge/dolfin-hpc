// Copyright (C) 2003 Johan Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_ARRAY_H
#define __DOLFIN_ARRAY_H

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>

#include <algorithm>
#include <iostream>
#include <vector>

namespace dolfin
{

template< typename T >
using Array = std::vector< T >;

template< typename T >
inline void free( Array< T * > & array )
{
  for ( T * element : array )
    if ( element != nullptr )
      delete element;

  array.clear();
}

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

} /* namespace dolfin */

#endif /* __DOLFIN_ARRAY_H */
