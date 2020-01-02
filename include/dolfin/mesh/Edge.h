// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_EDGE_H
#define __DOLFIN_EDGE_H

#include <dolfin/common/types.h>

#include "Point.h"
#include "MeshEntity.h"
#include "MeshEntityIterator.h"

namespace dolfin
{

/// An Edge is a MeshEntity of topological dimension 1.

class EdgeIterator;

class Edge : public MeshEntity
{

public:

  /// Create edge on given mesh
  Edge(Mesh& mesh, uint index) :
      MeshEntity(mesh, 1, index)
  {
  }

  /// Create edge from mesh entity
  Edge(MeshEntity& entity) :
      MeshEntity(entity.mesh(), 1, entity.index())
  {
  }

  /// Destructor
  ~Edge()
  {
  }

  /// Compute Euclidian length of edge
  real length() const;

  /// Compute coordinates of edge midpoint
  Point midpoint() const;

  //--- ITERATOR --------------------------------------------------------------

  typedef EdgeIterator iterator;

  struct shared : SharedIterator
  {
    shared(Mesh& M) : SharedIterator(M.topology().distdata()[1]) {}
    shared(MeshTopology& T) : SharedIterator(T.distdata()[1]) {}
  };

  struct ghost : GhostIterator
  {
    ghost(Mesh& M) : GhostIterator(M.topology().distdata()[1]) {}
    ghost(MeshTopology& T) : GhostIterator(T.distdata()[1]) {}
  };

  struct owned : OwnedIterator
  {
    owned(Mesh& M) : OwnedIterator(M.topology().distdata()[1]) {}
    owned(MeshTopology& T) : OwnedIterator(T.distdata()[1]) {}
  };

  //--- Entity relation -------------------------------------------------------

  typedef Vertex  lower_dimensional;
  typedef Face    higher_dimensional;

};

/// An EdgeIterator is a MeshEntityIterator of topological dimension 1.

class EdgeIterator : public MeshEntityIterator
{

public:

  EdgeIterator(Mesh& mesh) :
      MeshEntityIterator(mesh, 1)
  {
  }

  EdgeIterator(MeshEntity& entity) :
      MeshEntityIterator(entity, 1)
  {
  }

  inline Edge* operator->()
  {
    return static_cast<Edge*>(MeshEntityIterator::operator->());
  }

  inline Edge& operator*()
  {
    return *operator->();
  }

  inline Edge& operator[](uint i)
  {
    return static_cast<Edge&>(MeshEntityIterator::operator[](i));
  }

};

} /* namespace dolfin */

#endif /* __DOLFIN_EDGE_H */
