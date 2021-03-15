// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_EDGE_H
#define __DOLFIN_EDGE_H

#include <dolfin/common/GhostIterator.h>
#include <dolfin/common/OwnedIterator.h>
#include <dolfin/common/SharedIterator.h>
#include <dolfin/common/types.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/entities/MeshEntity.h>
#include <dolfin/mesh/entities/iterators/MeshEntityIterator.h>

namespace dolfin
{

class EdgeIterator;

//-----------------------------------------------------------------------------

/// An Edge is a MeshEntity of topological dimension 1.

class Edge : public MeshEntity
{

public:
  /// Create edge on given mesh
  Edge( Mesh & mesh, size_t index );

  /// Create edge from mesh entity
  Edge( MeshEntity & entity );

  /// Destructor
  ~Edge() = default;

  /// Compute Euclidian length of edge
  auto length() const -> real;

  /// Compute coordinates of edge midpoint
  auto midpoint() const -> Point;

  //--- ITERATOR --------------------------------------------------------------

  using iterator = EdgeIterator;

  struct shared : SharedIterator
  {
    shared( Mesh & M );
    shared( MeshTopology & T );
  };

  struct ghost : GhostIterator
  {
    ghost( Mesh & M );
    ghost( MeshTopology & T );
  };

  struct owned : OwnedIterator
  {
    owned( Mesh & M );
    owned( MeshTopology & T );
  };

  //--- Entity relation -------------------------------------------------------

  using lower_dimensional  = Vertex;
  using higher_dimensional = Face;
};

//-----------------------------------------------------------------------------

inline Edge::Edge( Mesh & mesh, size_t index )
  : MeshEntity( mesh, 1, index )
{
}

//-----------------------------------------------------------------------------

inline Edge::Edge( MeshEntity & entity )
  : MeshEntity( entity.mesh(), 1, entity.index() )
{
}

//-----------------------------------------------------------------------------

inline auto Edge::length() const -> real
{
  std::vector< size_t > const & vertices = entities( 0 );
  MeshGeometry const &          geom     = mesh_.geometry();

  dolfin_assert( not vertices.empty() );

  Point const & p0 = geom.point( vertices[0] );
  Point const & p1 = geom.point( vertices[1] );

  return p0.dist( p1 );
}

//-----------------------------------------------------------------------------

inline auto Edge::midpoint() const -> Point
{
  std::vector< size_t > const & vertices = entities( 0 );
  MeshGeometry const &          geom     = mesh_.geometry();

  dolfin_assert( not vertices.empty() );

  Point const & p0 = geom.point( vertices[0] );
  Point const & p1 = geom.point( vertices[1] );

  return 0.5 * ( p0 + p1 );
}

//-----------------------------------------------------------------------------

inline Edge::shared::shared( Mesh & M )
  : SharedIterator( M.topology().distdata()[1] )
{
}

//-----------------------------------------------------------------------------

inline Edge::shared::shared( MeshTopology & T )
  : SharedIterator( T.distdata()[1] )
{
}

//-----------------------------------------------------------------------------

inline Edge::ghost::ghost( Mesh & M )
  : GhostIterator( M.topology().distdata()[1] )
{
}

//-----------------------------------------------------------------------------

inline Edge::ghost::ghost( MeshTopology & T )
  : GhostIterator( T.distdata()[1] )
{
}

//-----------------------------------------------------------------------------

inline Edge::owned::owned( Mesh & M )
  : OwnedIterator( M.topology().distdata()[1] )
{
}

//-----------------------------------------------------------------------------

inline Edge::owned::owned( MeshTopology & T )
  : OwnedIterator( T.distdata()[1] )
{
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_EDGE_H */
