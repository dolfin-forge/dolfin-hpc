// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher 2017.
//
// First added:  2006-06-02
// Last changed: 2017-10-09

#ifndef __DOLFIN_FACE_H
#define __DOLFIN_FACE_H

#include "MeshEntity.h"
#include "MeshEntityIterator.h"

namespace dolfin
{

  /// A Face is a MeshEntity of topological dimension 2.

class FaceIterator;

class Face : public MeshEntity
{
public:

  /// Constructor
  Face(Mesh& mesh, uint index) : MeshEntity(mesh, 2, index) {}

  /// Destructor
  ~Face() {}

  /// Compute coordinates of face midpoint
  Point midpoint() const;

  //--- ITERATOR --------------------------------------------------------------

  typedef FaceIterator iterator;

  struct shared : SharedIterator
  {
    shared(MeshTopology& T) : SharedIterator(T.distdata()[2]) {}
  };

  struct ghost : GhostIterator
  {
    ghost(MeshTopology& T) : GhostIterator(T.distdata()[2]) {}
  };

  struct owned : OwnedIterator
  {
    owned(MeshTopology& T) : OwnedIterator(T.distdata()[2]) {}
  };


  //--- Entity relation -------------------------------------------------------

  typedef Edge  lower_dimensional;
  typedef Cell  higher_dimensional;

};

/// A FaceIterator is a MeshEntityIterator of topological dimension 2.

class FaceIterator : public MeshEntityIterator
{
public:
  
  FaceIterator(Mesh& mesh) : MeshEntityIterator(mesh, 2) {}
  FaceIterator(MeshEntity& entity) : MeshEntityIterator(entity, 2) {}

  inline Face& operator*() { return *operator->(); }
  inline Face* operator->() { return static_cast<Face*>(MeshEntityIterator::operator->()); }

};

} /* namespace dolfin */

#endif /* __DOLFIN_FACE_H */
