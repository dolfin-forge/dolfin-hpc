// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-06-21
// Last changed: 2008-05-28

#ifndef __BOUNDARY_MESH_H
#define __BOUNDARY_MESH_H

#include <dolfin/common/types.h>
#include "Mesh.h"
#include "MeshFunction.h"

namespace dolfin
{

/**
 *  @class  BoundaryMesh
 *
 *  @brief  A BoundaryMesh is a mesh over the boundary of some given mesh seen
 *          as a partition of a global mesh.
 */

class BoundaryMesh : public Mesh
{
public:

  enum Type { full, interior, exterior };

  /// Create boundary mesh from given mesh
  BoundaryMesh(Mesh& mesh, BoundaryMesh::Type type);

  /// Destructor
  ~BoundaryMesh();

  ///
  Mesh& mesh();

private:

  /// Global mesh hash
  Mesh& mesh_;
  std::string global_mesh_hash_;

};

}

#endif
