// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher 2017.
//
// First added:  2006-06-01
// Last changed: 2017-10-09

#ifndef __DOLFIN_VERTEX_H
#define __DOLFIN_VERTEX_H

#include "Point.h"
#include "MeshEntity.h"
#include "MeshEntityIterator.h"

namespace dolfin
{

/// A Vertex is a MeshEntity of topological dimension 0.

class VertexIterator;

class Vertex : public MeshEntity
{

public:

  /// Create vertex on given mesh
  Vertex(Mesh& mesh, uint index) :
      MeshEntity(mesh, 0, index)
  {
  }

  /// Create vertex from mesh entity
  Vertex(MeshEntity& entity) :
      MeshEntity(entity.mesh(), 0, entity.index())
  {
  }

  /// Destructor
  ~Vertex()
  {
  }

  /// Return array of vertex coordinates
  inline real* x()
  {
    return mesh_.geometry().x(index_);
  }

  /// Return array of vertex coordinates
  inline const real* x() const
  {
    return mesh_.geometry().x(index_);
  }

  /// Return vertex coordinates as a 3D point value
  inline Point point() const
  {
    return mesh_.geometry().point(index_);
  }

  //--- ITERATOR --------------------------------------------------------------

  typedef VertexIterator iterator;

  struct shared : SharedIterator
  {
    shared(MeshTopology& T) : SharedIterator(T.distdata()[0]) {}
  };

  struct ghost : GhostIterator
  {
    ghost(MeshTopology& T) : GhostIterator(T.distdata()[0]) {}
  };

  struct owned : OwnedIterator
  {
    owned(MeshTopology& T) : OwnedIterator(T.distdata()[0]) {}
  };

};

/// A VertexIterator is a MeshEntityIterator of topological dimension 0.

class VertexIterator : public MeshEntityIterator
{

public:

  VertexIterator(Mesh& mesh) :
      MeshEntityIterator(mesh, 0)
  {
  }

  VertexIterator(MeshEntity& entity) :
      MeshEntityIterator(entity, 0)
  {
  }

  inline Vertex* operator->()
  {
    return static_cast<Vertex*>(MeshEntityIterator::operator->());
  }
  
  inline Vertex& operator*()
  {
    return *operator->();
  }

  inline Vertex& operator[](uint i)
  {
    return static_cast<Vertex&>(MeshEntityIterator::operator[](i));
  }

};

} /* namespace dolfin */

#endif /* __DOLFIN_VERTEX_H */
