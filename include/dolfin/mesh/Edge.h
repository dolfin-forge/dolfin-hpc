// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_EDGE_H
#define __DOLFIN_EDGE_H

#include <dolfin/common/types.h>

#include "MeshEntity.h"
#include "MeshEntityIterator.h"
#include "Point.h"

#include <dolfin/common/GhostIterator.h>
#include <dolfin/common/OwnedIterator.h>
#include <dolfin/common/SharedIterator.h>

namespace dolfin
{

/// An Edge is a MeshEntity of topological dimension 1.

class EdgeIterator;

class Edge : public MeshEntity
{

public:
  /// Create edge on given mesh
  Edge( Mesh & mesh, uint index )
    : MeshEntity( mesh, 1, index )
  {
  }

  /// Create edge from mesh entity
  Edge( MeshEntity & entity )
    : MeshEntity( entity.mesh(), 1, entity.index() )
  {
  }

  /// Destructor
  ~Edge()
  {
  }

  /// Compute Euclidian length of edge
  real length() const;

  /// Compute coordinates of edge midpoint
  Point midpoint() const;

  //--- ITERATOR --------------------------------------------------------------

  typedef EdgeIterator iterator;

  struct shared : SharedIterator
  {
    shared( Mesh & M )
      : SharedIterator( M.topology().distdata()[1] )
    {
    }
    shared( MeshTopology & T )
      : SharedIterator( T.distdata()[1] )
    {
    }
  };

  struct ghost : GhostIterator
  {
    ghost( Mesh & M )
      : GhostIterator( M.topology().distdata()[1] )
    {
    }
    ghost( MeshTopology & T )
      : GhostIterator( T.distdata()[1] )
    {
    }
  };

  struct owned : OwnedIterator
  {
    owned( Mesh & M )
      : OwnedIterator( M.topology().distdata()[1] )
    {
    }
    owned( MeshTopology & T )
      : OwnedIterator( T.distdata()[1] )
    {
    }
  };

  //--- Entity relation -------------------------------------------------------

  typedef Vertex lower_dimensional;
  typedef Face   higher_dimensional;
};

inline real Edge::length() const
{
  Array< uint > const & vertices = entities( 0 );
  MeshGeometry const &  geom     = mesh_.geometry();

  dolfin_assert( not vertices.empty() );

  Point const & p0 = geom.point( vertices[0] );
  Point const & p1 = geom.point( vertices[1] );

  return p0.dist( p1 );
}
//-----------------------------------------------------------------------------
inline Point Edge::midpoint() const
{
  Array< uint > const & vertices = entities( 0 );
  MeshGeometry const &  geom     = mesh_.geometry();

  dolfin_assert( not vertices.empty() );

  Point const & p0 = geom.point( vertices[0] );
  Point const & p1 = geom.point( vertices[1] );

  return 0.5 * ( p0 + p1 );
}

} /* namespace dolfin */

#endif /* __DOLFIN_EDGE_H */
