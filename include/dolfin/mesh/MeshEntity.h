// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-11
// Last changed: 2006-10-23

#ifndef __DOLFIN_MESH_ENTITY_H
#define __DOLFIN_MESH_ENTITY_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

/**
 *
 *  @class  MeshEntity
 *
 *  @brief  A MeshEntity represents a mesh entity associated with a specific
 *          topological dimension of some mesh.
 *
 */

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

  //--- Topology --------------------------------------------------------------

  /// Return topological dimension
  uint dim() const;

  /// Return index of mesh entity
  uint index() const;

  /// Return number of incident mesh entities of given topological dimension
  uint num_entities(uint dim) const;

  /// Copy global indices of mesh entities to array
  void get_entities(uint dim, uint * indices) const;

  /// Copy global indices of mesh entities to array
  void get_entities(uint ** indices) const;

  /// Return array of indices for incident mesh entities of given topological
  /// dimension
  uint * entities(uint dim);

  /// Return array of indices for incident mesh entities of given topological
  /// dimension
  uint const * entities(uint dim) const;

  /// Check if given entity is incident
  bool incident(MeshEntity const& entity) const;

  /// Compute local index of given incident entity (-1 if not found)
  int index(MeshEntity const& entity) const;

  //--- Geometry --------------------------------------------------------------

  /*
   * TBD: add convenience function for accessing vertex coordinates
   *
   */

  //---------------------------------------------------------------------------

  /*
   * Convenience functions which can be called in serial and parallel.
   * For code portions within which the mesh is known to be distributed, calling
   * functions from the mesh distributed data is recommended.
   *
   */

  /// Return global index of mesh entity
  uint global_index() const;

  /// Copy global indices of mesh entities to array
  void get_global_entities(uint dim, uint * indices) const;

  /// Copy global indices of mesh entities to array
  void get_global_entities(uint ** indices) const;

  /// Return if the mesh entity is owned
  bool is_owned() const;

  /// Return if the mesh entity is shared
  bool is_shared() const;

  /// Return if the mesh entity is ghosted
  bool is_ghost() const;

  /// Return the owner of the mesh entity
  uint owner() const;

  /// Return if the mesh entity has all vertices shared
  bool has_all_vertices_shared() const;

  /// Return if the mesh entity is located on the global mesh boundary
  bool on_boundary() const;

  //---------------------------------------------------------------------------

  ///
  void disp() const;

  //--- CHECK ROUTINES --------------------------------------------------------

  /// Check
  void check() const;

protected:

  // Friends
  friend class MeshEntityIterator;

  // Mesh
  Mesh& mesh_;

  // Mesh topology
  MeshTopology& topology_;

  // Topological dimension
  uint const tdim_;

  // Geometric dimension
  uint const gdim_;

  // Pointer to mesh distributed data if applicable
  MeshDistributedData const * const distdata_;

  // Index of entity within topological dimension
  uint index_;

};

//--- INLINES -----------------------------------------------------------------

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
inline uint MeshEntity::num_entities(uint dim) const
{
  //dolfin_assert(topology_(tdim_, dim).is_initialized());
  //NOTE: New MeshTopology class auto-creates connectivity on demand.
  return topology_(tdim_, dim).degree(index_);
}

//-----------------------------------------------------------------------------
inline uint * MeshEntity::entities(uint dim)
{
  //dolfin_assert(topology_(tdim_, dim).is_initialized());
  //NOTE: New MeshTopology class auto-creates connectivity on demand.
  return topology_(tdim_, dim)(index_);
}

//-----------------------------------------------------------------------------
inline uint const * MeshEntity::entities(uint dim) const
{
  //dolfin_assert(topology_(tdim_, dim).is_initialized());
  //NOTE: New MeshTopology class auto-creates connectivity on demand.
  return topology_(tdim_, dim)(index_);
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_ENTITY_H */
