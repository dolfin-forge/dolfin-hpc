// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-03-01
// Last changed: 2007-03-01

#ifndef __DOLFIN_UFC_MESH_H
#define __DOLFIN_UFC_MESH_H

#include <ufc.h>

#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

/// This class is simple wrapper for a UFC mesh and provides
/// a layer between a DOLFIN mesh and a UFC mesh.

class UFCMesh : public ufc::mesh
{

public:

  /// Create UFC mesh from DOLFIN mesh
  UFCMesh(Mesh const& dolfin_mesh) :
      ufc::mesh(),
      mesh(&dolfin_mesh)
  {
    topological_dimension = mesh->topology_dimension();
    geometric_dimension = mesh->geometry().dim();
    num_entities = new uint[topological_dimension + 1];
    update();
  }

  /// Destructor
  ~UFCMesh()
  {
    delete[] num_entities;
    num_entities = NULL;
  }

  /// Create UFC mesh from DOLFIN mesh
  UFCMesh(UFCMesh const& other) :
      ufc::mesh(),
      mesh(other.mesh)
  {
    topological_dimension = mesh->topology_dimension();
    geometric_dimension = mesh->geometry().dim();
    num_entities = new uint[topological_dimension + 1];
    update();
  }

  ///
  Mesh const * const mesh;

  ///
  void update()
  {
    for (uint d = 0; d <= topological_dimension; ++d)
    {
      if (mesh->topology().entities_exist(d))
      {
        num_entities[d] = mesh->topology().global_size(d);
      }
    }
  }

};

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_MESH_H */
