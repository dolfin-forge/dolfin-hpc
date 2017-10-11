// Copyright (C) 2006 Johan Hoffman.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
//
// First added:  2006-11-01
// Last changed: 2010-09-13

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/LocalMeshRefinement.h>

#include <dolfin/main/MPI.h>
#include <dolfin/mesh/LoadBalancer.h>
#include <dolfin/mesh/RefinementManager.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace dolfin;

//-----------------------------------------------------------------------------
void LocalMeshRefinement::refineMeshByEdgeBisection(
    Mesh& mesh, MeshValues<bool, Cell>& cell_marker, bool refine_boundary, real tf,
    real tb, real ts, bool balance)
{
  begin("Refining simplicial mesh by edge bisection.");
  bool const is_distributed = mesh.is_distributed();

  // Start load balancer
  if (is_distributed)
  {
    // Generate cell - edge connectivity if not generated
    mesh.init(mesh.topology().dim(), 1);

    // Generate edge - vertex connectivity if not generated
    mesh.init(1, 0);

    if (balance)
    {
      begin("Load balancing");
      // Tune load balancer using machine specific parameters, if available
      if (tf > 0.0 && tb > 0.0 && ts > 0.0)
      {
        LoadBalancer::balance(mesh, cell_marker, tf, tb, ts);
      }
      else
      {
        LoadBalancer::balance(mesh, cell_marker);
      }
      end();
    }
  }

  // Get size of old mesh
  const uint num_vertices = mesh.size(0);
  const uint num_cells = mesh.size(mesh.topology().dim());

  // Check cell marker
  if (cell_marker.size() != num_cells) error("Wrong dimension of cell_marker");

  // Get cell type
  const CellType& cell_type = mesh.type();

  // Init new vertices and cells
  uint num_new_vertices = 0;
  uint num_new_cells = 0;

  // Create new mesh, refinement manager and open for editing
  Mesh refined_mesh;
  RefinementManager refman(mesh, refined_mesh);
  MeshEditor editor(refined_mesh, cell_type.cellType(), mesh.geometry().dim());

  // Initialize mappings
  Array<int> old2new_cell(mesh.num_cells());
  Array<int> old2new_vertex(mesh.size(0));

  // Initialise forbidden edges, default value is false
  MeshValues<bool, Edge> edge_forbidden(mesh);

  // If refinement of boundary is forbidden
  if (refine_boundary == false)
  {
    BoundaryMesh boundary(mesh, BoundaryMesh::exterior);
    // The boundary mesh could be empty in the parallel case
    if (boundary.size(0)) for (EdgeIterator e(boundary); !e.end(); ++e)
      edge_forbidden(e->index()) = true;
  }

  // Initialise forbidden cells
  MeshValues<bool, Cell> cell_forbidden(mesh);

  // Initialise forbidden cells on processor's boundary
  uint num_pcells = 0;
  uint num_pvertices = 0;

  if (is_distributed)
  {
    refman.mark_localboundary(cell_marker, num_pvertices, num_pcells);
  }

  // Initialise data for finding longest edge
  uint longest_edge_index = 0;
  real lmax, l;

  // Compute number of vertices and cells
  for (CellIterator c(mesh); !c.end(); ++c)
  {

    // Skip cell if marked as forbidden inside refinement manager
    if (refman.forbidden_cell(*c)) continue;

    if ((cell_marker(*c) == true) && (cell_forbidden(*c) == false))
    {
      // Find longest edge of cell c
      lmax = 0.0;
      for (EdgeIterator e(*c); !e.end(); ++e)
      {
        // Skip edges marked from propagation
        if (refman.forbidden_edge(*e)) continue;

        if (edge_forbidden(*e) == false)
        {
          l = e->length();
          if (lmax < l)
          {
            lmax = l;
            longest_edge_index = e->index();
          }
        }
      }

      Edge longest_edge(mesh, longest_edge_index);

      // If at least one edge should be bisected
      if (lmax > 0.0)
      {
        // Add new vertex
        num_new_vertices++;

        for (CellIterator cn(longest_edge); !cn.end(); ++cn)
        {
          if (cell_forbidden(*cn) == false)
          {
            // Count new cells
            num_new_cells++;
            // set markers of all cell neighbors of longest edge to false
            //if ( cn->index() != c->index() )
            cell_forbidden(cn->index()) = true;
            // set all the edges of cell neighbors to forbidden
            for (EdgeIterator en(*cn); !en.end(); ++en)
              edge_forbidden(en->index()) = true;
          }
        }
      }
    }
  }

  // Specify number of vertices and cells
  if (is_distributed)
  {
    editor.init_vertices(num_vertices + num_new_vertices + num_pvertices);
    editor.init_cells(num_cells + num_new_cells + num_pcells);
  }
  else
  {
    editor.init_vertices(num_vertices + num_new_vertices);
    editor.init_cells(num_cells + num_new_cells);
  }

  // Add old vertices
  uint current_vertex = 0;
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    if (is_distributed)
    {
      refined_mesh.distdata()[0].set_map(current_vertex, mesh.distdata()[0].get_global(v->index()));

      if (v->is_ghost())
      {
        refined_mesh.distdata()[0].set_ghost(current_vertex, mesh.distdata()[0].get_owner(v->index()));
      }
      else if (v->is_shared())
      {
        refined_mesh.distdata()[0].set_shared(current_vertex);
      }
    }

    editor.add_vertex(current_vertex++, v->point());
  }

  // Add old unrefined cells
  uint current_cell = 0;
  Array<uint> cell_vertices(cell_type.num_entities(0));
  for (CellIterator c(mesh); !c.end(); ++c)
  {

    // Skip unrefined cells which recives a propagated refinement
    if (refman.forbidden_cell(*c)) continue;

    //if ( (cell_marker(*c) == false) && (cell_forbidden(*c) == false) )
    if (cell_forbidden(*c) == false)
    {
      uint cv = 0;
      for (VertexIterator v(*c); !v.end(); ++v)
        cell_vertices[cv++] = v->index();
      editor.add_cell(current_cell++, &cell_vertices[0]);
    }

  }

  // Reset forbidden edges
  for (EdgeIterator e(mesh); !e.end(); ++e)
    edge_forbidden(e->index()) = false;

  // If refinement of boundary is forbidden
  if (refine_boundary == false)
  {
    BoundaryMesh boundary(mesh, BoundaryMesh::exterior);
    for (EdgeIterator e(boundary); !e.end(); ++e)
      edge_forbidden(e->index()) = true;
  }

  // Reset forbidden cells
  for (CellIterator c(mesh); !c.end(); ++c)
    cell_forbidden(c->index()) = false;

  // Add new vertices and cells.
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    if (is_distributed)
    {
      //      if( refman.forbidden_cell(*c) && !cell_forbidden(*c)) {
      if (refman.forbidden_cell(*c))
      {
        Edge e(mesh, refman.edge_refined(*c));
        if (edge_forbidden(e) == false)
        {
          edge_forbidden(e) = true;
          refman.add(e, current_vertex);
          editor.add_vertex(current_vertex++, e.midpoint());
          dolfin_assert(!cell_forbidden(*c));

          for (CellIterator cn(e); !cn.end(); ++cn)
          {
            dolfin_assert(!cell_forbidden(*cn));
            bisectEdgeOfSimplexCell(*cn, e, current_vertex, editor,
                                    current_cell);

            // Prevent any futher refinement on edge
            cell_marker(*cn) = false;
            cell_forbidden(cn->index()) = true;

          }
        }
        continue;
      }
    }

    if ((cell_marker(*c) == true) && (cell_forbidden(*c) == false))
    {
      // Find longest edge of cell c
      lmax = 0.0;
      for (EdgeIterator e(*c); !e.end(); ++e)
      {

        if (refman.forbidden_edge(*e)) continue;

        dolfin_assert(!refman.on_boundary(*e));

        if (edge_forbidden(*e) == false)
        {
          l = e->length();
          if (lmax < l)
          {
            lmax = l;
            longest_edge_index = e->index();
          }
        }
      }

      Edge longest_edge(mesh, longest_edge_index);

      // If at least one edge should be bisected
      if (lmax > 0.0)
      {
        dolfin_assert(!refman.on_boundary(longest_edge));
        refman.add(longest_edge, current_vertex);

        // Add new vertex
        editor.add_vertex(current_vertex++, longest_edge.midpoint());

        for (CellIterator cn(longest_edge); !cn.end(); ++cn)
        {
          dolfin_assert(!refman.forbidden_cell(*cn));

          // Add new cell
          bisectEdgeOfSimplexCell(*cn, longest_edge, current_vertex, editor,
                                  current_cell);

          // set markers of all cell neighbors of longest edge to false
          //if ( cn->index() != c->index() )
          cell_marker(cn->index()) = false;
          // set all edges of cell neighbors to forbidden
          for (EdgeIterator en(*cn); !en.end(); ++en)
            edge_forbidden(en->index()) = true;
        }
      }
    }
  }

  // Close edition and apply numbering of new entities
  editor.close();
  refman.apply();

  // Overwrite old mesh with refined mesh
  mesh = refined_mesh;
//  mesh.distdata().set_invalid_numbering();
//  mesh.renumber();

  end();
}
//-----------------------------------------------------------------------------
void LocalMeshRefinement::bisectEdgeOfSimplexCell(Cell& cell, Edge& edge,
                                                  uint& new_vertex,
                                                  MeshEditor& editor,
                                                  uint& current_cell)
{
  // Init cell vertices
  Array<uint> cell1_vertices(cell.num_entities(0));
  Array<uint> cell2_vertices(cell.num_entities(0));

  // Get edge vertices
  const uint* edge_vert = edge.entities(0);

  uint vc1 = 0;
  uint vc2 = 0;

  for (VertexIterator v(cell); !v.end(); ++v)
  {
    if ((v->index() != edge_vert[0]) && (v->index() != edge_vert[1]))
    {
      cell1_vertices[vc1++] = v->index();
      cell2_vertices[vc2++] = v->index();
    }
  }

  cell1_vertices[vc1++] = new_vertex - 1;
  cell2_vertices[vc2++] = new_vertex - 1;

  cell1_vertices[vc1++] = edge_vert[0];
  cell2_vertices[vc2++] = edge_vert[1];

  editor.add_cell(current_cell++, &cell1_vertices[0]);
  editor.add_cell(current_cell++, &cell2_vertices[0]);

}
//-----------------------------------------------------------------------------

