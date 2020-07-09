// Copyright (C) 2017. Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CLONABLE_H
#define __DOLFIN_CLONABLE_H

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace dolfin
{

template < class T >
struct Clonable
{
  /// Clone instance
  inline T * clone() const;

  /// Clone second argument into first argument if non-NULL
  static void clone( Clonable< T > *& p0, Clonable< T > const * p1 );
};

/// Clone if non-NULL, return NULL otherwise
template < class CloneT >
inline static CloneT * copyptr( CloneT const * p1 );

/// Clone if non-NULL, return NULL otherwise
template < class CloneT >
inline static CloneT * cloneptr( CloneT const * p1 );

//-----------------------------------------------------------------------------
template < class T >
inline T * Clonable< T >::clone() const
{
  return new T( static_cast< T const & >( *this ) );
}
//-----------------------------------------------------------------------------
template < class T >
void Clonable< T >::clone( Clonable< T > *& p0, Clonable< T > const * p1 )
{
  dolfin_assert( p0 == NULL );
  if ( p1 != NULL )
  {
    p0 = p1->clone();
  }
}
//-----------------------------------------------------------------------------
template < class CloneT >
inline static CloneT * copyptr( CloneT const * p1 )
{
  return ( p1 ? new CloneT( *p1 ) : NULL );
}
//-----------------------------------------------------------------------------
template < class CloneT >
inline static CloneT * cloneptr( CloneT const * p1 )
{
  return ( p1 ? p1->clone() : NULL );
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_CLONABLE_H */
