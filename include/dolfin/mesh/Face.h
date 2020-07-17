// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FACE_H
#define __DOLFIN_FACE_H

#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/MeshEntityIterator.h>

#include <dolfin/common/GhostIterator.h>
#include <dolfin/common/OwnedIterator.h>
#include <dolfin/common/SharedIterator.h>

namespace dolfin
{

/// A Face is a MeshEntity of topological dimension 2.

class FaceIterator;

class Face : public MeshEntity
{
public:
  /// Constructor
  Face( Mesh & mesh, uint index )
    : MeshEntity( mesh, 2, index )
  {
  }

  /// Destructor
  ~Face()
  {
  }

  /// Compute coordinates of face midpoint
  Point midpoint() const;

  //--- ITERATOR --------------------------------------------------------------

  using iterator = FaceIterator;

  struct shared : SharedIterator
  {
    shared( Mesh & M )
      : SharedIterator( M.topology().distdata()[2] )
    {
    }
    shared( MeshTopology & T )
      : SharedIterator( T.distdata()[2] )
    {
    }
  };

  struct ghost : GhostIterator
  {
    ghost( Mesh & M )
      : GhostIterator( M.topology().distdata()[2] )
    {
    }
    ghost( MeshTopology & T )
      : GhostIterator( T.distdata()[2] )
    {
    }
  };

  struct owned : OwnedIterator
  {
    owned( Mesh & M )
      : OwnedIterator( M.topology().distdata()[2] )
    {
    }
    owned( MeshTopology & T )
      : OwnedIterator( T.distdata()[2] )
    {
    }
  };

  //--- Entity relation -------------------------------------------------------

  using lower_dimensional  = Edge;
  using higher_dimensional = Cell;
};

inline Point Face::midpoint() const
{
  MeshGeometry const &  geometry     = this->mesh().geometry();
  Array< uint > const & vertices     = this->entities( 0 );
  uint const            num_vertices = this->num_entities( 0 );

  Point p;
  for ( uint v = 0; v < num_vertices; ++v )
  {
    p += geometry.point( vertices[v] );
  }
  p /= real( num_vertices );
  return p;
}

} /* namespace dolfin */

#endif /* __DOLFIN_FACE_H */
