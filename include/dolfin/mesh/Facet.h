// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FACET_H
#define __DOLFIN_FACET_H

#include <dolfin/common/DistributedData.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/MeshEntityIterator.h>

#include <dolfin/common/GhostIterator.h>
#include <dolfin/common/OwnedIterator.h>
#include <dolfin/common/SharedIterator.h>

namespace dolfin
{

/// A Facet is a MeshEntity of topological codimension 1.

class FacetIterator;

class Facet : public MeshEntity
{
public:
  /// Constructor
  Facet( Mesh & mesh, uint index )
    : MeshEntity( mesh, mesh.type().facet_dim(), index )
  {
  }

  /// Destructor
  ~Facet() = default;

  /// Compute coordinates of facet midpoint
  Point midpoint() const;

  //--- ITERATOR --------------------------------------------------------------

  using iterator = FacetIterator;

  struct shared : SharedIterator
  {
    shared( Mesh & M )
      : SharedIterator( M.topology().distdata()[M.type().facet_dim()] )
    {
    }
    shared( MeshTopology & T )
      : SharedIterator( T.distdata()[T.type().facet_dim()] )
    {
    }
  };

  struct ghost : GhostIterator
  {
    ghost( Mesh & M )
      : GhostIterator( M.topology().distdata()[M.type().facet_dim()] )
    {
    }
    ghost( MeshTopology & T )
      : GhostIterator( T.distdata()[T.type().facet_dim()] )
    {
    }
  };

  struct owned : OwnedIterator
  {
    owned( Mesh & M )
      : OwnedIterator( M.topology().distdata()[M.type().facet_dim()] )
    {
    }
    owned( MeshTopology & T )
      : OwnedIterator( T.distdata()[T.type().facet_dim()] )
    {
    }
  };
};

//-----------------------------------------------------------------------------
inline Point Facet::midpoint() const
{
  MeshGeometry const &  geometry     = this->mesh().geometry();
  Array< uint > const & vertices     = this->entities( 0 );
  uint const            num_vertices = this->num_entities( 0 );

  Point p;
  for ( uint v = 0; v < num_vertices; ++v )
  {
    p += geometry.point( vertices[v] );
  }
  p /= static_cast< real >( num_vertices );
  return p;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_FACET_H */
