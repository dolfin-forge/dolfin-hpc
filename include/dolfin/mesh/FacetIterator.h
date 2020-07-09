// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FACET_ITERATOR_H
#define __DOLFIN_FACET_ITERATOR_H

#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/MeshEntityIterator.h>

namespace dolfin
{

/// A FacetIterator is a MeshEntityIterator of topological codimension 1.

class FacetIterator : public MeshEntityIterator
{
public:
  FacetIterator( Mesh & mesh )
    : MeshEntityIterator( mesh, mesh.type().facet_dim() )
  {
  }
  FacetIterator( MeshEntity & entity )
    : MeshEntityIterator( entity, entity.mesh().type().facet_dim() )
  {
  }

  inline Facet * operator->();

  inline Facet & operator*();

  inline Facet & operator[]( uint i );
};

inline Facet * FacetIterator::operator->()
{
  return static_cast< Facet * >( MeshEntityIterator::operator->() );
}

inline Facet & FacetIterator::operator*()
{
  return *operator->();
}

inline Facet & FacetIterator::operator[]( uint i )
{
  return static_cast< Facet & >( MeshEntityIterator::operator[]( i ) );
}

} /* namespace dolfin */

#endif /* __DOLFIN_FACET_ITERATOR_H */
