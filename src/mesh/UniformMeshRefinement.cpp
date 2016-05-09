// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
// Modified by Stefanie Strunk, 2013.
//
// First added:  2006-06-08
// Last changed: 2008-07-07

#include <dolfin/mesh/UniformMeshRefinement.h>

#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/RefinementManager.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
void UniformMeshRefinement::refine(Mesh& mesh)
{
  message(1, "Refining %s mesh uniformly.", mesh.type().str().c_str());

  // Create new mesh, refinement manager and open for editing
  uint const tdim = mesh.topology().dim();
  Mesh refined_mesh;
  MeshEditor editor(refined_mesh, mesh.type(), mesh.geometry().dim());
  RefinementManager refman(mesh, refined_mesh);
  RefinementPattern const& pattern = refman.pattern();

  // Compute number of vertices based on the pattern
  uint num_refined_vertices = mesh.topology().size(0);
  for (uint i = 1; i <= tdim; ++i)
  {
    dolfin_assert(mesh.topology().size(i) > 0);
    num_refined_vertices += pattern.num_refined_vertices(i)
        * mesh.topology().size(i);
  }
  editor.init_vertices(num_refined_vertices);

  // Refinement pattern provides the number of new cells
  uint num_refined_cells = pattern.num_refined_cells() * mesh.num_cells();
  editor.init_cells(num_refined_cells);

  // Current vertex index
  uint current_vertex = 0;

  // Add old vertices
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    refman.add(*v, current_vertex);
    editor.add_vertex(current_vertex++, v->point());
  }

  // Add edge-based vertices
  if (tdim > 1 && pattern.refinement_needs_entities(1))
  {
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      refman.add(*e, current_vertex);
      editor.add_vertex(current_vertex++, e->midpoint());
    }
  }

  // Add face-based vertices
  if (tdim > 2 && pattern.refinement_needs_entities(2))
  {
    for (FaceIterator f(mesh); !f.end(); ++f)
    {
      refman.add(*f, current_vertex);
      editor.add_vertex(current_vertex++, f->midpoint());
    }
  }

  // Add cell-based vertices
  if (pattern.refinement_needs_entities(tdim))
  {
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      editor.add_vertex(current_vertex++, c->midpoint());
    }
  }

  // Add cells
  uint current_cell = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    pattern.refine_cell(*c, editor, current_cell);
  }

  // Apply numbering of new entities and close edition
  refman.apply();
  editor.close();

  // Overwrite old mesh with refined mesh
  mesh = refined_mesh;
  mesh.topology().renumber();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

