// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_VERTEX_ITERATOR_H
#define __DOLFIN_VERTEX_ITERATOR_H

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/MeshEntity.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/MeshEntityIterator.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

/**
 *  @class  VertexIterator
 *
 *  @brief  A VertexIterator is a MeshEntityIterator of topological dimension 0.
 *
 */

class VertexIterator : public MeshEntityIterator
{
public:
  VertexIterator( Mesh & mesh );
  VertexIterator( MeshEntity & entity );

  inline auto operator->() -> Vertex *;
  inline auto operator*() -> Vertex &;
  inline auto operator[]( size_t i ) -> Vertex &;
};

//-----------------------------------------------------------------------------

inline VertexIterator::VertexIterator( Mesh & mesh )
  : MeshEntityIterator( mesh, 0 )
{
}

//-----------------------------------------------------------------------------

inline VertexIterator::VertexIterator( MeshEntity & entity )
  : MeshEntityIterator( entity, 0 )
{
}

inline auto VertexIterator::operator->() -> Vertex *
{
  return static_cast< Vertex * >( MeshEntityIterator::operator->() );
}

//-----------------------------------------------------------------------------

inline auto VertexIterator::operator*() -> Vertex &
{
  return *operator->();
}

//-----------------------------------------------------------------------------

inline auto VertexIterator::operator[]( size_t i ) -> Vertex &
{
  return static_cast< Vertex & >( MeshEntityIterator::operator[]( i ) );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_VERTEX_ITERATOR_H */
