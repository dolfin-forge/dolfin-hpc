// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

// This class was rewritten to support hypercube cells and fixes several issues:
// - invalid computation of the boundary mesh with missing or added cells.
// - segmentation faults due to non-robust handling of zero-size boundaries.

// The reimplementation creates boundaries of distributed meshes which are
// themselves distributed, allowing nested calls.

#ifndef __DOLFIN_BOUNDARY_MESH_H
#define __DOLFIN_BOUNDARY_MESH_H

#include <dolfin/common/Array.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshDependent.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

class SubDomain;

/**
 *  @class  BoundaryMesh
 *
 *  @brief  A BoundaryMesh is a mesh over the boundary of some given mesh seen
 *          as a partition of a global mesh.
 */

class BoundaryMesh : public Mesh, public MeshDependent
{

public:

  enum Type { full, interior, exterior };

  /// Create boundary mesh from given mesh
  BoundaryMesh(Mesh& mesh, BoundaryMesh::Type type);

  /// Create boundary mesh from given boundary mesh
  explicit BoundaryMesh(BoundaryMesh& mesh, BoundaryMesh::Type type);

  /// Create boundary mesh from given mesh
  BoundaryMesh(Mesh& mesh, SubDomain const& subdomain, BoundaryMesh::Type type);

  /// Create boundary mesh from subdomain of given boundary mesh
  BoundaryMesh(BoundaryMesh& boundary, SubDomain const& subdomain, bool inside);

  /// Copy constructor
  BoundaryMesh(BoundaryMesh const& other);

  /// Destructor
  ~BoundaryMesh() override;

  /// Return facet index in the mesh associated with the boundary cell
  uint facet_index(Cell const& boundary_cell) const;

  /// Return facet index in the mesh associated with the boundary cell
  uint facet_index(uint boundary_cell_index) const;

  /// Return vertex index in the mesh associated with the boundary vertex
  uint vertex_index(Vertex const& boundary_vertex) const;

  /// Return vertex index in the mesh associated with the boundary vertex
  uint vertex_index(uint boundary_vertex_index) const;

  /// Return type
  Type boundary_type() const;

  /// Return whether it is the boundary of a boundary mesh
  bool is_boundary_of_boundary() const;

  /// Return if the topology should be reordered
  bool reordering() const override { return false; }

private:

  ///
  void init(Mesh& mesh, BoundaryMesh::Type type);

  ///
  void compute(Mesh& mesh, bool exterior, bool interior);

  ///
  Type type_;

  ///
  bool const boundary_of_boundary_;

  ///
  Array<uint> cell_map_;
  Array<uint> vertex_map_;

  ///
  SubDomain const * subdomain_;

};

//-----------------------------------------------------------------------------
inline uint BoundaryMesh::facet_index(Cell const& boundary_cell) const
{
  dolfin_assert(&boundary_cell.mesh() == this);
  return cell_map_[boundary_cell.index()];
}
//-----------------------------------------------------------------------------
inline uint BoundaryMesh::facet_index(uint boundary_cell_index) const
{
  dolfin_assert(boundary_cell_index < cell_map_.size());
  return cell_map_[boundary_cell_index];
}
//-----------------------------------------------------------------------------
inline uint BoundaryMesh::vertex_index(Vertex const& boundary_vertex) const
{
  dolfin_assert(&boundary_vertex.mesh() == this);
  return vertex_map_[boundary_vertex.index()];
}
//-----------------------------------------------------------------------------
inline uint BoundaryMesh::vertex_index(uint boundary_vertex_index) const
{
  dolfin_assert(boundary_vertex_index < vertex_map_.size());
  return vertex_map_[boundary_vertex_index];
}
//-----------------------------------------------------------------------------
inline BoundaryMesh::Type BoundaryMesh::boundary_type() const
{
  return type_;
}
//-----------------------------------------------------------------------------
inline bool BoundaryMesh::is_boundary_of_boundary() const
{
  return boundary_of_boundary_;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_BOUNDARY_MESH_H */
