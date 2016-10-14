// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2016.
//
// This class was rewritten to support hypercube cells and fixes several issues:
// - invalid computation of the boundary mesh with missing or added cells.
// - segmentation faults due to non-robust handling of zero-size boundaries.
//
// The reimplementation creates boundaries of distributed meshes which are
// themselves distributed, allowing nested calls.
//
// First added:  2006-06-21
// Last changed: 2008-05-28

#ifndef __DOLFIN_BOUNDARY_MESH_H
#define __DOLFIN_BOUNDARY_MESH_H

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshDependent.h>

#include <dolfin/common/Array.h>

namespace dolfin
{

class Cell;
class SubDomain;
class Vertex;

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

  /// Create boundary mesh from given mesh
  BoundaryMesh(Mesh& mesh, SubDomain const& subdomain, BoundaryMesh::Type type);

  /// Create boundary mesh from subdomain of given boundary mesh
  BoundaryMesh(BoundaryMesh& boundary, SubDomain const& subdomain, bool inside);

  /// Destructor
  ~BoundaryMesh();

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

private:

  ///
  void compute(Mesh& mesh, bool exterior, bool interior);

  ///
  Type type_;

  ///
  Array<uint> cell_map_;
  Array<uint> vertex_map_;

  ///
  SubDomain const * subdomain_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_BOUNDARY_MESH_H */
