// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_VERTEX_H
#define __DOLFIN_VERTEX_H

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/Point.h>

#include <dolfin/common/GhostIterator.h>
#include <dolfin/common/OwnedIterator.h>
#include <dolfin/common/SharedIterator.h>

namespace dolfin
{

/**
 *  @class  Vertex
 *
 *  @brief  A Vertex is a MeshEntity of topological dimension 0.
 *
 */

class VertexIterator;

class Vertex : public MeshEntity
{

public:
  /// Create vertex on given mesh
  Vertex( Mesh & mesh, uint index )
    : MeshEntity( mesh, 0, index )
  {
  }

  /// Create vertex from mesh entity
  Vertex( MeshEntity & entity )
    : MeshEntity( entity.mesh(), 0, entity.index() )
  {
  }

  /// Destructor
  ~Vertex() = default;

  /// Return array of vertex coordinates
  inline auto x() -> real *;

  /// Return array of vertex coordinates
  inline auto x() const -> const real *;

  /// Return vertex coordinates as a 3D point value
  inline auto point() const -> Point;

  /// Return vertex coordinates as a 3D point value
  inline auto midpoint() const -> Point;

  //--- ITERATOR --------------------------------------------------------------

  using iterator = VertexIterator;

  struct shared : SharedIterator
  {
    shared( Mesh & M )
      : SharedIterator( M.topology().distdata()[0] )
    {
    }
    shared( MeshTopology & T )
      : SharedIterator( T.distdata()[0] )
    {
    }
  };

  struct ghost : GhostIterator
  {
    ghost( Mesh & M )
      : GhostIterator( M.topology().distdata()[0] )
    {
    }
    ghost( MeshTopology & T )
      : GhostIterator( T.distdata()[0] )
    {
    }
  };

  struct owned : OwnedIterator
  {
    owned( Mesh & M )
      : OwnedIterator( M.topology().distdata()[0] )
    {
    }
    owned( MeshTopology & T )
      : OwnedIterator( T.distdata()[0] )
    {
    }
  };

  //--- Entity relation -------------------------------------------------------

  using lower_dimensional  = Vertex;
  using higher_dimensional = Edge;
};

//-----------------------------------------------------------------------------
inline auto Vertex::x() -> real *
{
  return mesh_.geometry().x( index_ );
}

//-----------------------------------------------------------------------------
inline auto Vertex::x() const -> const real *
{
  return mesh_.geometry().x( index_ );
}

//-----------------------------------------------------------------------------
inline auto Vertex::point() const -> Point
{
  return mesh_.geometry().point( index_ );
}

//-----------------------------------------------------------------------------
inline auto Vertex::midpoint() const -> Point
{
  return mesh_.geometry().point( index_ );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_VERTEX_H */
