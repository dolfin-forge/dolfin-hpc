// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
//
// First added:  2006-06-08
// Last changed: 2007-05-24

#include <dolfin/math/dolfin_math.h>
#include <dolfin/log/dolfin_log.h>
#include "Mesh.h"
#include "MeshTopology.h"
#include "MeshGeometry.h"
#include "MeshConnectivity.h"
#include "MeshEditor.h"
#include "Vertex.h"
#include "Edge.h"
#include "Cell.h"
#include "UniformMeshRefinement.h"

#include "RefinementManager.h"
#include "MeshFunction.h"
#include <dolfin/main/MPI.h>
#include <mpi.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
void UniformMeshRefinement::refine(Mesh& mesh)
{
  // Only know how to refine simplicial meshes
  refineSimplex(mesh);
}
//-----------------------------------------------------------------------------
void UniformMeshRefinement::refineSimplex(Mesh& mesh)
{
  message(1, "Refining simplicial mesh uniformly.");
  
  // Generate cell - edge connectivity if not generated
  mesh.init(mesh.topology().dim(), 1);
  
  // Generate edge - vertex connectivity if not generated
  mesh.init(1, 0);

  // Get cell type
  const CellType& cell_type = mesh.type();

  // Create a new refinement manager  FIXME:remove for the serial case
  RefinementManager refman(mesh);

  // Create new mesh and open for editing
  Mesh refined_mesh;
  MeshEditor editor;
  editor.open(refined_mesh, cell_type.cellType(),
	      mesh.topology().dim(), mesh.geometry().dim());
  
  // Get size of mesh
  const uint num_vertices = mesh.size(0);
  const uint num_edges = mesh.size(1);
  const uint num_cells = mesh.size(mesh.topology().dim());

  // Specify number of vertices and cells
  editor.initVertices(num_vertices + num_edges);
  editor.initCells(ipow(2, mesh.topology().dim())*num_cells);

  uint* edge_vert;
  Array<uint> shared_edge;
    
  uint vertex = 0;
  if(MPI::numProcesses() > 1){
    // Add old vertices
    for (VertexIterator v(mesh); !v.end(); ++v) {
      refined_mesh.distdata().set_map(vertex, mesh.distdata().get_global(*v), 0);
      
      if( mesh.distdata().is_ghost(v->index()) ) {
	refined_mesh.distdata().set_ghost(vertex);
	refined_mesh.distdata().set_ghost_owner(vertex,
						mesh.distdata().get_owner(*v));
      }
      else if(mesh.distdata().is_shared(v->index()))
	refined_mesh.distdata().set_shared(vertex);

      editor.addVertex(vertex++, v->point());
    }
    
    for (EdgeIterator e(mesh); !e.end(); ++e) {
      edge_vert = e->entities(0);
      // If the edge is shared and lies between processes
      // process new vertex inside refinement manager
      if( refman.on_boundary(*e) ){
	
	// Add the new vertex inside the refinement manager
	refman.addVertex(edge_vert, vertex, refined_mesh);
	
	// Buffer edge information for mapping phase
	shared_edge.push_back(edge_vert[0]);
	shared_edge.push_back(edge_vert[1]);
	shared_edge.push_back(vertex);
      } 
      else
        refman.addVertex(vertex, refined_mesh);

      editor.addVertex(vertex++, e->midpoint());
    }

  }
  else {
    
    // Add old vertices
    for (VertexIterator v(mesh); !v.end(); ++v)
      editor.addVertex(vertex++, v->point());
    
    // Add new vertices
    for (EdgeIterator e(mesh); !e.end(); ++e)
      editor.addVertex(vertex++, e->midpoint());
    
  }

  // Add cells
  uint current_cell = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
    cell_type.refineCell(*c, editor, current_cell);

  editor.close();
  
  // Map global numbers to unassigned shared vertices
  if(MPI::numProcesses() > 1)
    refman.map_new_vertices(shared_edge, mesh, refined_mesh);

  // Overwrite old mesh with refined mesh
  mesh = refined_mesh;

  if(MPI::numProcesses() >1) {
    mesh.distdata().invalid_numbering();
    // FIXME, fix map_new_vertices such that all datastructures are working
    MeshFunction<dolfin::uint> part;
    part.init(mesh, mesh.topology().dim());
    part = dolfin::MPI::processNumber();
    mesh.distribute(part);
    mesh.renumber();
    message("%d", mesh.distdata().global_numVertices());

    uint tmp =  mesh.numVertices() - mesh.distdata().num_ghost();
    // MPI aliasing 
    uint num_glb;  
    MPI_Allreduce(&tmp, &num_glb, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    
    mesh.distdata().set_global_numVertices(num_glb);
  }

}
//-----------------------------------------------------------------------------
