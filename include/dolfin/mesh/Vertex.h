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
  ~Vertex()
  {
  }

  /// Return array of vertex coordinates
  inline real * x();

  /// Return array of vertex coordinates
  inline const real * x() const;

  /// Return vertex coordinates as a 3D point value
  inline Point point() const;

  /// Return vertex coordinates as a 3D point value
  inline Point midpoint() const;

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
inline real * Vertex::x()
{
  return mesh_.geometry().x( index_ );
}

//-----------------------------------------------------------------------------
inline const real * Vertex::x() const
{
  return mesh_.geometry().x( index_ );
}

//-----------------------------------------------------------------------------
inline Point Vertex::point() const
{
  return mesh_.geometry().point( index_ );
}

//-----------------------------------------------------------------------------
inline Point Vertex::midpoint() const
{
  return mesh_.geometry().point( index_ );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_VERTEX_H */
