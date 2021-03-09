// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FACET_H
#define __DOLFIN_FACET_H

#include <dolfin/common/DistributedData.h>
#include <dolfin/common/GhostIterator.h>
#include <dolfin/common/OwnedIterator.h>
#include <dolfin/common/SharedIterator.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/MeshEntity.h>
#include <dolfin/mesh/entities/iterators/MeshEntityIterator.h>

namespace dolfin
{

class FacetIterator;

//-----------------------------------------------------------------------------

/// A Facet is a MeshEntity of topological codimension 1.

class Facet : public MeshEntity
{
public:
  /// Constructor
  Facet( Mesh & mesh, size_t index );

  /// Destructor
  ~Facet() = default;

  /// Compute coordinates of facet midpoint
  auto midpoint() const -> Point;

  //--- ITERATOR --------------------------------------------------------------

  using iterator = FacetIterator;

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
};

//-----------------------------------------------------------------------------

inline Facet::Facet( Mesh & mesh, size_t index )
  : MeshEntity( mesh, mesh.type().facet_dim(), index )
{
}

//-----------------------------------------------------------------------------

inline auto Facet::midpoint() const -> Point
{
  MeshGeometry const &          geometry     = this->mesh().geometry();
  std::vector< size_t > const & vertices     = this->entities( 0 );
  size_t const                  num_vertices = this->num_entities( 0 );

  Point p;
  for ( size_t v = 0; v < num_vertices; ++v )
  {
    p += geometry.point( vertices[v] );
  }
  p /= static_cast< real >( num_vertices );
  return p;
}

//-----------------------------------------------------------------------------

inline Facet::shared::shared( Mesh & M )
  : SharedIterator( M.topology().distdata()[M.type().facet_dim()] )
{
}

//-----------------------------------------------------------------------------

inline Facet::shared::shared( MeshTopology & T )
  : SharedIterator( T.distdata()[T.type().facet_dim()] )
{
}

//-----------------------------------------------------------------------------

inline Facet::ghost::ghost( Mesh & M )
  : GhostIterator( M.topology().distdata()[M.type().facet_dim()] )
{
}

//-----------------------------------------------------------------------------

inline Facet::ghost::ghost( MeshTopology & T )
  : GhostIterator( T.distdata()[T.type().facet_dim()] )
{
}

//-----------------------------------------------------------------------------

inline Facet::owned::owned( Mesh & M )
  : OwnedIterator( M.topology().distdata()[M.type().facet_dim()] )
{
}

//-----------------------------------------------------------------------------

inline Facet::owned::owned( MeshTopology & T )
  : OwnedIterator( T.distdata()[T.type().facet_dim()] )
{
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_FACET_H */
