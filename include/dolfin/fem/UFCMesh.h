// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-03-01
// Last changed: 2007-03-01

#ifndef __DOLFIN_UFC_MESH_H
#define __DOLFIN_UFC_MESH_H

#include <ufc.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

/// This class is simple wrapper for a UFC mesh and provides
/// a layer between a DOLFIN mesh and a UFC mesh.

class UFCMesh : public ufc::mesh
{
public:

  /// Create emtpy UFC mesh
  UFCMesh() :
      ufc::mesh(),
      mesh(NULL)
  {
  }

  /// Create UFC mesh from DOLFIN mesh
  UFCMesh(Mesh& dolfin_mesh) :
      ufc::mesh(),
      mesh(&dolfin_mesh)
  {
    init(dolfin_mesh);
  }

  /// Destructor
  ~UFCMesh()
  {
    clear();
  }

  Mesh const * mesh;

  /// Initialize UFC mesh data
  void init(Mesh const& mesh);

  /// Display info
  void disp() const;

private:

  // Clear UFC cell data
  void clear();

};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline void UFCMesh::init(Mesh const& mesh)
{
  clear();

  // Update pointer to current mesh
  this->mesh = &mesh;

  // Set topological dimension
  uint const tdim = mesh.topology().dim();
  topological_dimension = tdim;

  // Set geometric dimension
  geometric_dimension = mesh.geometry().dim();

  // Set number of entities for each topological dimension
  num_entities = new uint[tdim + 1];
  for(uint d = 0; d < tdim + 1; ++d)
  {
    num_entities[d] = mesh.topology().num_global(d);
  }
}

//-----------------------------------------------------------------------------
inline void UFCMesh::clear()
{
  mesh = NULL;
  topological_dimension = 0;
  geometric_dimension = 0;
  delete[] num_entities;
  num_entities = NULL;
}

//-----------------------------------------------------------------------------
inline void UFCMesh::disp() const
{
  cout << "UFCMesh" << endl;
  cout << "-------" << endl;

  // Begin indentation
  begin("");

  // Display UFC dofmap information
  cout << "ufc::mesh info" << endl;
  cout << "--------------" << endl;
  begin("");
  cout << "Topological dimension : " << topological_dimension << endl;
  cout << "Geometric dimension   : " << geometric_dimension << endl;
  for (uint d = 0; d <= topological_dimension; ++d)
  {
    cout << "Number of entities of dim(" << d << ") : " << num_entities[d]
         << endl;
  }
  cout << endl;
  end();

  // End indentation
  end();
}

}

#endif
