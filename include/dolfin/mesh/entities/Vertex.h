// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_VERTEX_H
#define __DOLFIN_VERTEX_H

#include <dolfin/common/GhostIterator.h>
#include <dolfin/common/OwnedIterator.h>
#include <dolfin/common/SharedIterator.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/entities/MeshEntity.h>

namespace dolfin
{

class VertexIterator;

//-----------------------------------------------------------------------------

/**
 *  @class  Vertex
 *
 *  @brief  A Vertex is a MeshEntity of topological dimension 0.
 *
 */

class Vertex : public MeshEntity
{

public:
  /// Create vertex on given mesh
  Vertex( Mesh & mesh, size_t index );

  /// Create vertex from mesh entity
  Vertex( MeshEntity & entity );

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
  using higher_dimensional = Edge;
};

//-----------------------------------------------------------------------------

inline Vertex::Vertex( Mesh & mesh, size_t index )
  : MeshEntity( mesh, 0, index )
{
}

//-----------------------------------------------------------------------------

/// Create vertex from mesh entity
inline Vertex::Vertex( MeshEntity & entity )
  : MeshEntity( entity.mesh(), 0, entity.index() )
{
}

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

inline Vertex::shared::shared( Mesh & M )
  : SharedIterator( M.topology().distdata()[0] )
{
}

//-----------------------------------------------------------------------------

inline Vertex::shared::shared( MeshTopology & T )
  : SharedIterator( T.distdata()[0] )
{
}

//-----------------------------------------------------------------------------

inline Vertex::ghost::ghost( Mesh & M )
  : GhostIterator( M.topology().distdata()[0] )
{
}

//-----------------------------------------------------------------------------

inline Vertex::ghost::ghost( MeshTopology & T )
  : GhostIterator( T.distdata()[0] )
{
}

//-----------------------------------------------------------------------------

inline Vertex::owned::owned( Mesh & M )
  : OwnedIterator( M.topology().distdata()[0] )
{
}

//-----------------------------------------------------------------------------

inline Vertex::owned::owned( MeshTopology & T )
  : OwnedIterator( T.distdata()[0] )
{
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_VERTEX_H */
