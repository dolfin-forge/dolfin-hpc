// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-06-21
// Last changed: 2008-05-28

#ifndef __BOUNDARY_MESH_H
#define __BOUNDARY_MESH_H

#include <dolfin/common/types.h>

#include "Mesh.h"
#include "MeshDependent.h"

namespace dolfin
{

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

  /// Destructor
  ~BoundaryMesh();

  /// Return facet index in the mesh associated with the boundary cell
  uint facet_index(Cell const& boundary_cell);

  /// Return vertex index in the mesh associated with the boundary vertex
  uint vertex_index(Vertex const& boundary_vertex);

private:


};

}

#endif
