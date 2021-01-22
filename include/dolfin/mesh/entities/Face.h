// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FACE_H
#define __DOLFIN_FACE_H

#include <dolfin/common/GhostIterator.h>
#include <dolfin/common/OwnedIterator.h>
#include <dolfin/common/SharedIterator.h>
#include <dolfin/mesh/entities/MeshEntity.h>
#include <dolfin/mesh/entities/iterators/MeshEntityIterator.h>

namespace dolfin
{

class FaceIterator;

//-----------------------------------------------------------------------------

/// A Face is a MeshEntity of topological dimension 2.

class Face : public MeshEntity
{
public:
  /// Constructor
  Face( Mesh & mesh, size_t index );

  /// Destructor
  ~Face() = default;

  /// Compute coordinates of face midpoint
  auto midpoint() const -> Point;

  //--- ITERATOR --------------------------------------------------------------

  using iterator = FaceIterator;

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

  using lower_dimensional  = Edge;
  using higher_dimensional = Cell;
};

//-----------------------------------------------------------------------------

inline Face::Face( Mesh & mesh, size_t index )
  : MeshEntity( mesh, 2, index )
{
}

//-----------------------------------------------------------------------------

inline auto Face::midpoint() const -> Point
{
  MeshGeometry const &          geometry     = this->mesh().geometry();
  std::vector< size_t > const & vertices     = this->entities( 0 );
  size_t const                  num_vertices = this->num_entities( 0 );

  Point p;
  for ( size_t v = 0; v < num_vertices; ++v )
  {
    p += geometry.point( vertices[v] );
  }
  p /= real( num_vertices );
  return p;
}

//-----------------------------------------------------------------------------

inline Face::shared::shared( Mesh & M )
  : SharedIterator( M.topology().distdata()[2] )
{
}

//-----------------------------------------------------------------------------

inline Face::shared::shared( MeshTopology & T )
  : SharedIterator( T.distdata()[2] )
{
}

//-----------------------------------------------------------------------------

inline Face::ghost::ghost( Mesh & M )
  : GhostIterator( M.topology().distdata()[2] )
{
}

//-----------------------------------------------------------------------------

inline Face::ghost::ghost( MeshTopology & T )
  : GhostIterator( T.distdata()[2] )
{
}

//-----------------------------------------------------------------------------

inline Face::owned::owned( Mesh & M )
  : OwnedIterator( M.topology().distdata()[2] )
{
}

//-----------------------------------------------------------------------------

inline Face::owned::owned( MeshTopology & T )
  : OwnedIterator( T.distdata()[2] )
{
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_FACE_H */
