// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2012-2016.
//
// This class was initially cured in 2012 from a terrible disease "infinite
// loop on zero-size connectivity".
// The absence of robustness of the original class is fixed by "on-demand"
// computation of connectivities in the mesh topology.
//
// First added:  2006-05-09
// Last changed: 2007-05-02

#ifndef __DOLFIN_MESH_ENTITY_ITERATOR_H
#define __DOLFIN_MESH_ENTITY_ITERATOR_H

#include <dolfin/common/types.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEntity.h>

namespace dolfin
{

/// MeshEntityIterator provides a common iterator for mesh entities
/// over meshes, boundaries and incidence relations. The basic use
/// is illustrated below.
///
/// The following example shows how to iterate over all mesh entities
/// of a mesh of topological dimension dim:
///
///     for (MeshEntityIterator e(mesh, dim); !e.end(); ++e)
///     {
///       e->foo();
///     }
///
/// The following example shows how to iterate over mesh entities of
/// topological dimension dim connected (incident) to some mesh entity f:
///
///     for (MeshEntityIterator e(f, dim); !e.end(); ++e)
///     {
///       e->foo();
///     }
///
/// In addition to the general iterator, a set of specific named iterators
/// are provided for entities of type Vertex, Edge, Face, Facet and Cell.
/// These iterators are defined along with their respective classes.

class MeshEntityIterator
{
public:

  /// Create iterator for mesh entities over given topological dimension
  MeshEntityIterator(Mesh& mesh, uint dim);

  /// Create iterator for entities of given dimension connected to given entity
  MeshEntityIterator(MeshEntity& entity, uint dim);

  /// Destructor
  virtual ~MeshEntityIterator();

  /// Step to next mesh entity (prefix increment)
  MeshEntityIterator& operator++();

  /// Return current position
  uint pos() const;

  /// Check if iterator has reached the end
  bool end() const;

  /// Dereference operator
  MeshEntity& operator*();

  /// Member access operator
  MeshEntity* operator->();

private:

  /// Copy constructor is private to disallow usage. If it were public (or not
  /// declared and thus a default version available) it would allow code like
  ///
  /// for (CellIterator c0(mesh); !c0.end(); ++c0)
  ///   for (CellIterator c1(c0); !c1.end(); ++c1)
  ///      ...
  ///
  /// c1 looks to be an iterator over the entities around c0 when it is in
  /// fact a copy of c0.
  MeshEntityIterator(MeshEntityIterator& other);

  // Mesh entity
  MeshEntity entity_;

  // Current position
  uint pos_;

  // End position
  uint const pos_end_;

  // Mapping from pos to index (if any)
  uint const * index_;
  
};

//--- INLINES -----------------------------------------------------------------

inline MeshEntityIterator& MeshEntityIterator::operator++()
{
  ++pos_;
  return *this;
}

//-----------------------------------------------------------------------------
inline uint MeshEntityIterator::pos() const
{
  return pos_;
}

//-----------------------------------------------------------------------------
inline bool MeshEntityIterator::end() const
{
  return pos_ >= pos_end_;
}

//-----------------------------------------------------------------------------
inline MeshEntity& MeshEntityIterator::operator*()
{
  return *operator->();
}

//-----------------------------------------------------------------------------
inline MeshEntity* MeshEntityIterator::operator->()
{
  entity_.index_ = (index_ ? index_[pos_] : pos_);
  return &entity_;
}

//-----------------------------------------------------------------------------

}

#endif
