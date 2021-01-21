// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_VERTEX_ITERATOR_H
#define __DOLFIN_VERTEX_ITERATOR_H

#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/MeshEntityIterator.h>

namespace dolfin
{

/**
 *  @class  VertexIterator
 *
 *  @brief  A VertexIterator is a MeshEntityIterator of topological dimension 0.
 *
 */

class VertexIterator : public MeshEntityIterator
{
public:
  VertexIterator( Mesh & mesh )
    : MeshEntityIterator( mesh, 0 )
  {
  }

  VertexIterator( MeshEntity & entity )
    : MeshEntityIterator( entity, 0 )
  {
  }

inline auto operator->() -> Vertex *
{
  return static_cast< Vertex * >( MeshEntityIterator::operator->() );
}

inline auto operator*() -> Vertex &
{
  return *operator->();
}

inline auto operator[]( uint i ) -> Vertex &
{
  return static_cast< Vertex & >( MeshEntityIterator::operator[]( i ) );
}

};



} /* namespace dolfin */

#endif /* __DOLFIN_VERTEX_ITERATOR_H */
