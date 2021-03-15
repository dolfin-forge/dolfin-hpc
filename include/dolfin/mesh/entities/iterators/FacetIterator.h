// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FACET_ITERATOR_H
#define __DOLFIN_FACET_ITERATOR_H

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/MeshEntity.h>
#include <dolfin/mesh/entities/iterators/MeshEntityIterator.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

/// A FacetIterator is a MeshEntityIterator of topological codimension 1.

class FacetIterator : public MeshEntityIterator
{
public:
  FacetIterator( Mesh & mesh );
  FacetIterator( MeshEntity & entity );

  inline auto operator->() -> Facet *;
  inline auto operator*() -> Facet &;
  inline auto operator[]( size_t i ) -> Facet &;
};

//-----------------------------------------------------------------------------

inline FacetIterator::FacetIterator( Mesh & mesh )
  : MeshEntityIterator( mesh, mesh.type().facet_dim() )
{
}

//-----------------------------------------------------------------------------

inline FacetIterator::FacetIterator( MeshEntity & entity )
  : MeshEntityIterator( entity, entity.mesh().type().facet_dim() )
{
}

inline auto FacetIterator::operator->() -> Facet *
{
  return static_cast< Facet * >( MeshEntityIterator::operator->() );
}

//-----------------------------------------------------------------------------

inline auto FacetIterator::operator*() -> Facet &
{
  return *operator->();
}

//-----------------------------------------------------------------------------

inline auto FacetIterator::operator[]( size_t i ) -> Facet &
{
  return static_cast< Facet & >( MeshEntityIterator::operator[]( i ) );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_FACET_ITERATOR_H */
