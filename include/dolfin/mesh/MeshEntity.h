// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-11
// Last changed: 2006-10-23

#ifndef __MESH_ENTITY_H
#define __MESH_ENTITY_H

#include <dolfin/common/types.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

/// A MeshEntity represents a mesh entity associated with
/// a specific topological dimension of some mesh.

class MeshEntity
{
public:

  /// Constructor
  MeshEntity(Mesh& mesh, uint dim, uint index);

  /// Destructor
  ~MeshEntity();

  /// Return mesh associated with mesh entity
  Mesh& mesh();

  /// Return mesh associated with mesh entity
  Mesh const& mesh() const;

  /// Return topological dimension
  uint dim() const;

  /// Return index of mesh entity
  uint index() const;

  /// Return number of incident mesh entities of given topological dimension
  uint numEntities(uint dim) const;

  /// Return array of indices for incident mesh entities of given topological dimension
  uint * entities(uint dim);

  /// Return array of indices for incident mesh entities of given topological dimension
  uint const * entities(uint dim) const;

  /// Check if given entity is incident
  bool incident(MeshEntity const& entity) const;

  /// Compute local index of given incident entity (error if not found)
  uint index(MeshEntity const& entity) const;

  /// Return if the mesh entity is shared
  bool is_shared() const;

  /// Return if the mesh entity is ghosted
  bool is_ghost() const;

  /// Output
  friend LogStream& operator<<(LogStream& stream, MeshEntity const& entity);

protected:

  // Friends
  friend class MeshEntityIterator;

  // The mesh
  Mesh& mesh_;

  // Topological dimension
  uint const tdim_;

  // Geometric dimension
  uint const gdim_;

  // Index of entity within topological dimension
  uint index_;

};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline Mesh& MeshEntity::mesh()
{
  return mesh_;
}

//-----------------------------------------------------------------------------
inline Mesh const& MeshEntity::mesh() const
{
  return mesh_;
}

//-----------------------------------------------------------------------------
inline uint MeshEntity::dim() const
{
  return tdim_;
}

//-----------------------------------------------------------------------------
inline uint MeshEntity::index() const
{
  return index_;
}

//-----------------------------------------------------------------------------
inline uint MeshEntity::numEntities(uint dim) const
{
  return mesh_.topology()(tdim_, dim).size(index_);
}

//-----------------------------------------------------------------------------
inline uint * MeshEntity::entities(uint dim)
{
  return mesh_.topology()(tdim_, dim)(index_);
}

//-----------------------------------------------------------------------------
inline uint const * MeshEntity::entities(uint dim) const
{
  return mesh_.topology()(tdim_, dim)(index_);
}

}

#endif
