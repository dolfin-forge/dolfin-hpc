// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
// Modified by Stefanie Strunk, 2013.
//
// First added:  2006-06-08
// Last changed: 2008-07-07

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
#include <dolfin/mesh/UniformMeshRefinement.h>

#include <dolfin/mesh/RefinementManager.h>
#include <dolfin/main/MPI.h>

#ifdef HAVE_LIBGEOM
#include <Geometry.h>
#endif

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------
void UniformMeshRefinement::refine(Mesh& mesh)
{
  switch (mesh.type().cellType())
    {
    case CellType::interval:
    case CellType::triangle:
    case CellType::tetrahedron:
      refineSimplex(mesh);
      break;
    case CellType::quadrilateral:
    case CellType::hexahedron:
      refineHypercube(mesh);
      break;
    default:
      error("Uniform mesh refinement not implemented for given mesh type.");
      break;
    }
}
//-----------------------------------------------------------------------------
void UniformMeshRefinement::refineSimplex(Mesh& mesh)
{
  message(1, "Refining simplicial mesh uniformly.");

  uint const tdim = mesh.topology().dim();

  // Create new mesh, refinement manager and open for editing
  Mesh refined_mesh;
  CellType const& cell_type = mesh.type();
  MeshEditor editor(refined_mesh, cell_type.cellType(), mesh.geometry().dim());
  RefinementManager refman(mesh, refined_mesh);
  RefinementPattern const& pattern = refman.pattern();

  // Specify number of vertices and cells
  editor.init_vertices(mesh.size(0) + mesh.size(1));
  editor.init_cells(pattern.num_refined_cells() * mesh.size(tdim));

  // Current vertex index
  uint vertex = 0;

  // Add old vertices
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    refman.add(*v, vertex);
    editor.add_vertex(vertex++, v->point());
  }

  // Add edge-based vertices
  for (EdgeIterator e(mesh); !e.end(); ++e)
  {
    refman.add(*e, vertex);
    editor.add_vertex(vertex++, e->midpoint());
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
  //mesh.renumber();
}
//-----------------------------------------------------------------------------
void UniformMeshRefinement::refineHypercube(Mesh& mesh)
{
  message(1, "Refining n-cube mesh uniformly.");

  uint const tdim = mesh.topology().dim();

  // Create new mesh, refinement manager and open for editing
  Mesh refined_mesh;
  CellType const& cell_type = mesh.type();
  MeshEditor editor(refined_mesh, cell_type.cellType(), mesh.geometry().dim());
  RefinementManager refman(mesh, refined_mesh);
  RefinementPattern const& pattern = refman.pattern();

  // Refinement pattern creates one new vertex per entity of each dimension
  uint num_refined_vertices = mesh.topology().size(0);
  for (uint i = 1; i <= tdim; ++i)
  {
    dolfin_assert(mesh.topology().size(i) > 0);
    num_refined_vertices += mesh.topology().size(i);
  }
  editor.init_vertices(num_refined_vertices);

  // Refinement pattern provides the number of new cells
  uint num_refined_cells = cell_type.num_refined_cells() * mesh.num_cells();
  editor.init_cells(num_refined_cells);

  // Current vertex index
  uint vertex = 0;

  // Add vertices for each topological dimension
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    refman.add(*v, vertex);
    editor.add_vertex(vertex++, v->point());
  }

  // Add edge-based vertices
  if (tdim > 1)
  {
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      refman.add(*e, vertex);
      editor.add_vertex(vertex++, e->midpoint());
    }
  }

  // Add face-based vertices
  if (tdim > 2)
  {
    for (FaceIterator f(mesh); !f.end(); ++f)
    {
      refman.add(*f, vertex);
      editor.add_vertex(vertex++, f->midpoint());
    }
  }

  // Add cell-based vertices and cells
  uint cell = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    editor.add_vertex(vertex++, c->midpoint());
    pattern.refine_cell(*c, editor, cell);
  }

  // Apply numbering of new entities and close edition
  refman.apply();
  editor.close();

  // Overwrite old mesh with refined mesh
  mesh = refined_mesh;
  //mesh.renumber();
}
//-----------------------------------------------------------------------------

}

