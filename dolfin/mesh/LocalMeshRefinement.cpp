// Copyright (C) 2006 Johan Hoffman.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-11-01
// Last changed: 2007-01-16

#include <dolfin/math/dolfin_math.h>
#include <dolfin/log/dolfin_log.h>
#include "Mesh.h"
#include "MeshTopology.h"
#include "MeshGeometry.h"
#include "MeshConnectivity.h"
#include "MeshEditor.h"
#include "MeshFunction.h"
#include "Vertex.h"
#include "Edge.h"
#include "Cell.h"
#include "BoundaryMesh.h"
#include "LocalMeshRefinement.h"

#include <dolfin/main/MPI.h>
#include "LoadBalancer.h"
#include "RefinementManager.h"

#include <mpi.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
void LocalMeshRefinement::refineMeshByEdgeBisection(Mesh& mesh, 
                                                    MeshFunction<bool>& cell_marker, 
						    bool refine_boundary, 
						    real tf, real tb, real ts)
{
  begin("Refining simplicial mesh by edge bisection.");

  // Start Loadbalanacer
  if(MPI::numProcesses() > 1) {
    // Generate cell - edge connectivity if not generated
    mesh.init(mesh.topology().dim(), 1);
    
    // Generate edge - vertex connectivity if not generated
    mesh.init(1, 0);
    
    begin("Load balancing");
    // Tune loadbalancer using machine specific parameters if available
    if( tf > 0.0 && tb > 0.0 && ts > 0.0)
      LoadBalancer::balance(mesh, cell_marker, tf, tb, ts);
    else
      LoadBalancer::balance(mesh, cell_marker);
    end();
  }

  // Get size of old mesh
  const uint num_vertices = mesh.size(0);
  const uint num_cells = mesh.size(mesh.topology().dim());
  
  // Check cell marker 
  if ( cell_marker.size() != num_cells ) error("Wrong dimension of cell_marker");
  
  // Generate cell - edge connectivity if not generated
  mesh.init(mesh.topology().dim(), 1);
  
  // Generate edge - vertex connectivity if not generated
  mesh.init(1, 0);
  
  // Get cell type
  const CellType& cell_type = mesh.type();

  // Create a new refinement manager
  RefinementManager refman(mesh); 
  Array<uint> shared_edge;
  
  // Init new vertices and cells
  uint num_new_vertices = 0;
  uint num_new_cells = 0;
  
  // Compute number of vertices and cells 
  /*
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    if((cell_marker.get(*c) == true))
    {
      //cout << "marked cell: " << endl;
      //cout << c->midpoint() << endl;
      //cout << c->index() << endl;
    }
  }
  */

  // Create new mesh and open for editing
  Mesh refined_mesh;
  MeshEditor editor;
  editor.open(refined_mesh, cell_type.cellType(),
	      mesh.topology().dim(), mesh.geometry().dim());
  
  // Initialize mappings
  Array<int> old2new_cell(mesh.numCells());
  Array<int> old2new_vertex(mesh.numVertices());

  // Initialise forbidden edges 
  MeshFunction<bool> edge_forbidden(mesh);  
  edge_forbidden.init(1);
  for (EdgeIterator e(mesh); !e.end(); ++e)
    edge_forbidden.set(e->index(),false);
  
  // If refinement of boundary is forbidden 
  if ( refine_boundary == false )
  {
    BoundaryMesh boundary(mesh);
    // The boundary mesh could be empty in the parallel case
    if( boundary.size(0) )
      for (EdgeIterator e(boundary); !e.end(); ++e)
	edge_forbidden.set(e->index(),true);
  }

  // Initialise forbidden cells 
  MeshFunction<bool> cell_forbidden(mesh);  
  cell_forbidden.init(mesh.topology().dim());
  for (CellIterator c(mesh); !c.end(); ++c)
    cell_forbidden.set(c->index(),false);
  
  // Initialise forbidden cells on process boundary
  uint num_propagated, num_pcells, num_nv;
  if(MPI::numProcesses() > 1) 
    refman.mark_localboundary(mesh, cell_marker, num_nv, num_propagated, num_pcells);

  // Initialise data for finding longest edge   
  uint longest_edge_index = 0;
  real lmax, l;

  // Compute number of vertices and cells 
  for (CellIterator c(mesh); !c.end(); ++c)
  {

    // Skip cell if marked as forbidden inside refinement manager
    if( refman.forbidden_cell(*c) )
      continue;
    
    if ( (cell_marker.get(*c) == true) && (cell_forbidden.get(*c) == false) )
    {
//       cout << "marked cell: " << endl;
//       cout << c->midpoint() << endl;
//       cout << c->index() << endl;


      // Find longest edge of cell c
      lmax = 0.0;
      for (EdgeIterator e(*c); !e.end(); ++e)
      {
	// Skip edges marked from propagation
	if( refman.forbidden_edge(*e) )
	  continue;
	
	if ( edge_forbidden.get(*e) == false )
	{
	  l = e->length();
	  if ( lmax < l )
	  {
	    lmax = l;
	    longest_edge_index = e->index(); 
	  }
	}
      }

      Edge longest_edge(mesh,longest_edge_index);

      // If at least one edge should be bisected
      if ( lmax > 0.0 )
      {
	// Add new vertex 
	num_new_vertices++;

	for (CellIterator cn(longest_edge); !cn.end(); ++cn)
	{
          if ( cell_forbidden.get(*cn) == false )
          {
            // Count new cells
            num_new_cells++;
            // set markers of all cell neighbors of longest edge to false 
            //if ( cn->index() != c->index() ) 
            cell_forbidden.set(cn->index(),true);
            // set all the edges of cell neighbors to forbidden
            for (EdgeIterator en(*cn); !en.end(); ++en)
              edge_forbidden.set(en->index(),true);
          }
	}
      }
    }
  }
  
  uint nv = 0;
  // Specify number of vertices and cells
  if(MPI::numProcesses() > 1) {
    editor.initVertices(num_vertices + num_new_vertices + num_propagated + num_nv);
    editor.initCells(num_cells + num_new_cells + num_pcells);
    nv = num_vertices + num_new_vertices + num_propagated + num_nv;
  }
  else {
    editor.initVertices(num_vertices + num_new_vertices);
    editor.initCells(num_cells + num_new_cells);
  }

  //cout << "Number of cells in old mesh: " << num_cells << "; to add: " << num_new_cells << endl;
  //cout << "Number of vertices in old mesh: " << num_vertices << "; to add: " << num_new_vertices << endl;
  
  // Add old vertices
  uint current_vertex = 0;
  for (VertexIterator v(mesh); !v.end(); ++v) {
    if(MPI::numProcesses() > 1) {
      refined_mesh.distdata().set_map(current_vertex, 
				      mesh.distdata().get_global(*v), 0);
      
      if(mesh.distdata().is_ghost(v->index())) {
	refined_mesh.distdata().set_ghost(current_vertex);
	refined_mesh.distdata().set_ghost_owner(current_vertex, mesh.distdata().get_owner(v->index()));
      }
      else if(mesh.distdata().is_shared(v->index()))
	refined_mesh.distdata().set_shared(current_vertex);
    }
      //    else
      editor.addVertex(current_vertex++, v->point());

  }

  // Add old unrefined cells 
  uint current_cell = 0;
  Array<uint> cell_vertices(cell_type.numEntities(0));
  for (CellIterator c(mesh); !c.end(); ++c)
  {
  
    // Skip unrefined cells which recives a propagated refinement
    if(refman.forbidden_cell(*c))
      continue;

    //if ( (cell_marker.get(*c) == false) && (cell_forbidden.get(*c) == false) )
    if ( cell_forbidden.get(*c) == false )
    {
      uint cv = 0;
      for (VertexIterator v(*c); !v.end(); ++v)
        cell_vertices[cv++] = v->index(); 
      editor.addCell(current_cell++, cell_vertices);
    }

  }
  
  // Reset forbidden edges 
  for (EdgeIterator e(mesh); !e.end(); ++e)
    edge_forbidden.set(e->index(),false);

  // If refinement of boundary is forbidden 
  if ( refine_boundary == false )
  {
    BoundaryMesh boundary(mesh);
    for (EdgeIterator e(boundary); !e.end(); ++e)
      edge_forbidden.set(e->index(),true);
  }

  // Reset forbidden cells 
  for (CellIterator c(mesh); !c.end(); ++c)
    cell_forbidden.set(c->index(),false);
  
  shared_edge.clear();

  // Add new vertices and cells. 
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    if(MPI::numProcesses() > 1) {
      if( refman.forbidden_cell(*c) && !cell_forbidden.get(*c)) {      
	Edge e(mesh, refman.edge_refined(*c));
	if(edge_forbidden.get(e) == false) {
	  edge_forbidden.set(e, true);
	  uint *edge_vert = e.entities(0);
	  shared_edge.push_back(edge_vert[0]);
	  shared_edge.push_back(edge_vert[1]);
	  shared_edge.push_back(current_vertex);

	  refman.addVertex(edge_vert, current_vertex, refined_mesh);      
	  editor.addVertex(current_vertex++, e.midpoint());
	  dolfin_assert( !cell_forbidden.get(*c) );

	  
	  for(CellIterator cn(e); !cn.end(); ++cn) {
	    dolfin_assert( !cell_forbidden.get(*cn) );
	    bisectEdgeOfSimplexCell(*cn, e, current_vertex, editor, current_cell);
	    
	    // Prevent any futher refinement on edge
	    cell_marker.set(*cn, false);
	    cell_forbidden.set(cn->index(), true);
	    
	  }
	}
	continue;
      }
    }

    if ( (cell_marker.get(*c) == true) && (cell_forbidden.get(*c) == false) )
    {
      // Find longest edge of cell c
      lmax = 0.0;
      for (EdgeIterator e(*c); !e.end(); ++e)
      {

	if( refman.forbidden_edge(*e) )
	  continue;
	if( refman.on_boundary(*e))
	  error("Missing forbidden marking on shared_edge");

	if ( edge_forbidden.get(*e) == false )
	{
	  l = e->length();
	  if ( lmax < l )
	  {
	    lmax = l;
	    longest_edge_index = e->index(); 
	  }
	}
      }
     
      Edge longest_edge(mesh,longest_edge_index);

      // If at least one edge should be bisected
      if ( lmax > 0.0 )
      {
	dolfin_assert( !refman.on_boundary(longest_edge) );
	refman.addVertex(current_vertex, refined_mesh);

	// Add new vertex
	editor.addVertex(current_vertex++, longest_edge.midpoint());

	//cout << "adding new vertex: " << endl;
	//cout << longest_edge.midpoint() << endl;

	for (CellIterator cn(longest_edge); !cn.end(); ++cn)
	{
	  // Add new cell
	  bisectEdgeOfSimplexCell(*cn, longest_edge, current_vertex, editor, current_cell);

	  // set markers of all cell neighbors of longest edge to false 
	  //if ( cn->index() != c->index() ) 
          cell_marker.set(cn->index(),false);
	  // set all edges of cell neighbors to forbidden
	  for (EdgeIterator en(*cn); !en.end(); ++en)
	    edge_forbidden.set(en->index(),true);
	}
      }
    }
  }

  // Assign global numbers to shared vertices
  if(MPI::numProcesses() > 1) 
    refman.map_new_vertices(shared_edge, mesh, refined_mesh);

  // Overwrite old mesh with refined mesh
  editor.close();
  mesh = refined_mesh;

  if(MPI::numProcesses() >1) {
    // FIXME, fix map_new_vertices such that all datastructures are working
    MeshFunction<dolfin::uint> part;
    part.init(mesh, mesh.topology().dim());
    part = dolfin::MPI::processNumber();
    mesh.distribute(part);
    mesh.renumber();
  }
  end();
}
//-----------------------------------------------------------------------------
void LocalMeshRefinement::bisectEdgeOfSimplexCell(Cell& cell, 
                                                  Edge& edge, 
                                                  uint& new_vertex, 
                                                  MeshEditor& editor, 
                                                  uint& current_cell) 
{
  // Init cell vertices 
  Array<uint> cell1_vertices(cell.numEntities(0));
  Array<uint> cell2_vertices(cell.numEntities(0));

  // Get edge vertices 
  const uint* edge_vert = edge.entities(0);

  uint vc1 = 0;
  uint vc2 = 0;

  for (VertexIterator v(cell); !v.end(); ++v)
  {
    if ( (v->index() != edge_vert[0]) && (v->index() != edge_vert[1]) )
    {
      cell1_vertices[vc1++] = v->index();
      cell2_vertices[vc2++] = v->index();
    }
  }

  cell1_vertices[vc1++] = new_vertex - 1; 
  cell2_vertices[vc2++] = new_vertex - 1; 
  
  cell1_vertices[vc1++] = edge_vert[0]; 
  cell2_vertices[vc2++] = edge_vert[1]; 

  editor.addCell(current_cell++, cell1_vertices);
  editor.addCell(current_cell++, cell2_vertices);

}
//-----------------------------------------------------------------------------

