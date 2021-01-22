// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_EDGE_ITERATOR_H
#define __DOLFIN_EDGE_ITERATOR_H

#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/entities/Edge.h>
#include <dolfin/mesh/entities/MeshEntity.h>
#include <dolfin/mesh/entities/iterators/MeshEntityIterator.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

/// An EdgeIterator is a MeshEntityIterator of topological dimension 1.

class EdgeIterator : public MeshEntityIterator
{

public:
  EdgeIterator( Mesh & mesh );
  EdgeIterator( MeshEntity & entity );

  inline auto operator->() -> Edge *;
  inline auto operator*() -> Edge &;
  inline auto operator[]( uint i ) -> Edge &;
};

//-----------------------------------------------------------------------------

inline EdgeIterator::EdgeIterator( Mesh & mesh )
  : MeshEntityIterator( mesh, 1 )
{
}

//-----------------------------------------------------------------------------

inline EdgeIterator::EdgeIterator( MeshEntity & entity )
  : MeshEntityIterator( entity, 1 )
{
}

//-----------------------------------------------------------------------------

inline auto EdgeIterator::operator->() -> Edge *
{
  return static_cast< Edge * >( MeshEntityIterator::operator->() );
}

//-----------------------------------------------------------------------------

inline auto EdgeIterator::operator*() -> Edge &
{
  return *operator->();
}

//-----------------------------------------------------------------------------

inline auto EdgeIterator::operator[]( uint i ) -> Edge &
{
  return static_cast< Edge & >( MeshEntityIterator::operator[]( i ) );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_EDGE_ITERATOR_H */
