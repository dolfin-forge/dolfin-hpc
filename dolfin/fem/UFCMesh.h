// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-03-01
// Last changed: 2007-03-01

#ifndef __UFC_MESH_H
#define __UFC_MESH_H

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
    
    /// Create empty UFC mesh
    UFCMesh() : ufc::mesh() {}

    /// Create UFC mesh from DOLFIN mesh
    UFCMesh(Mesh& mesh) : ufc::mesh()
    {
      init(mesh);
    }

    /// Destructor
    ~UFCMesh()
    {
      clear();
    }

    /// Initialize UFC cell data
    void init(Mesh& mesh)
    {
      // Clear old data
      clear();

      // Set topological dimension
      topological_dimension = mesh.topology().dim();
      
      // Set geometric dimension
      geometric_dimension = mesh.geometry().dim();

      // Set number of entities for each topological dimension
      num_entities = new uint[mesh.topology().dim() + 1];
      for (uint d = 0; d <= mesh.topology().dim(); d++)
	if( d == 0 && MPI::numProcesses() > 1)
	  num_entities[0] = mesh.distdata().global_numVertices();
	else if( d == 1 && MPI::numProcesses() > 1)
	  num_entities[1] = mesh.distdata().global_numEdges();
	else if( d == 2 && MPI::numProcesses() > 1)
	  if( mesh.topology().dim() > 2)
	    num_entities[2] = mesh.distdata().global_numFaces();
	  else
	    num_entities[2] = mesh.distdata().global_numCells();
	else
	  num_entities[d] = mesh.size(d);
    }

    // Clear UFC cell data
    void clear()
    {
      topological_dimension = 0;
      geometric_dimension = 0;

      if ( num_entities )
        delete [] num_entities;
      num_entities = 0;
    }

  };

}

#endif
