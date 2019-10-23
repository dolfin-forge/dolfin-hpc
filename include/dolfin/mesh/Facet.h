// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FACET_H
#define __DOLFIN_FACET_H

#include "MeshEntity.h"
#include "MeshEntityIterator.h"

namespace dolfin
{

/// A Facet is a MeshEntity of topological codimension 1.

class FacetIterator;

class Facet : public MeshEntity
{
public:

  /// Constructor
  Facet(Mesh& mesh, uint index) : MeshEntity(mesh, mesh.type().facet_dim(), index) {}

  /// Destructor
  ~Facet() {}

  /// Compute coordinates of facet midpoint
  Point midpoint() const;

  //--- ITERATOR --------------------------------------------------------------

  typedef FacetIterator iterator;

  struct shared : SharedIterator
  {
    shared(Mesh& M) : SharedIterator(M.topology().distdata()[M.type().facet_dim()]) {}
    shared(MeshTopology& T) : SharedIterator(T.distdata()[T.type().facet_dim()]) {}
  };

  struct ghost : GhostIterator
  {
    ghost(Mesh& M) : GhostIterator(M.topology().distdata()[M.type().facet_dim()]) {}
    ghost(MeshTopology& T) : GhostIterator(T.distdata()[T.type().facet_dim()]) {}
  };

  struct owned : OwnedIterator
  {
    owned(Mesh& M) : OwnedIterator(M.topology().distdata()[M.type().facet_dim()]) {}
    owned(MeshTopology& T) : OwnedIterator(T.distdata()[T.type().facet_dim()]) {}
  };

};

/// A FacetIterator is a MeshEntityIterator of topological codimension 1.

class FacetIterator : public MeshEntityIterator
{
public:
  
  FacetIterator(Mesh& mesh) : MeshEntityIterator(mesh, mesh.type().facet_dim()) {}
  FacetIterator(MeshEntity& entity) : MeshEntityIterator(entity, entity.mesh().type().facet_dim()) {}

  inline Facet* operator->() { return static_cast<Facet*>(MeshEntityIterator::operator->()); }
  inline Facet& operator*() { return *operator->(); }
  inline Facet& operator[](uint i) { return static_cast<Facet&>(MeshEntityIterator::operator[](i)); }

};

} /* namespace dolfin */

#endif /* __DOLFIN_FACET_H */
