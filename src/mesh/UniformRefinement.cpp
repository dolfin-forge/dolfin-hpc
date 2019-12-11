// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/UniformRefinement.h>

#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Face.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Connectivity.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
template<class E>
void add_refined_vertices(MeshEditor& editor, Mesh& mesh)
{
  Mesh& refined_mesh = editor.mesh();
  uint const tdim = refined_mesh.topology_dimension();
  uint const edim = entity_dimension<E>(mesh);
  if (tdim > edim && mesh.type().num_refined_vertices(edim))
  {
    uint const voffset = editor.current_vertex();
    for (typename E::iterator e(mesh); !e.end(); ++e)
    {
      editor.add_vertex(voffset + e->index(), e->midpoint());
    }
    if (mesh.is_distributed())
    {
      uint goffset = 0;
      for (uint i = 0; i < edim; ++i)
      {
        goffset += mesh.global_size(i) * mesh.type().num_refined_vertices(i);
      }
      DistributedData& dist = refined_mesh.distdata()[0];
      for (typename E::iterator e(mesh); !e.end(); ++e)
      {
        dist.set_map(voffset + e->index(), goffset + e->global_index());
      }
      for (typename E::shared e(mesh); e.valid(); ++e)
      {
        dist.setall_shared_adj(voffset + e.index(), e.adj());
      }
      for (typename E::ghost e(mesh); e.valid(); ++e)
      {
        dist.set_ghost(voffset + e.index(), e.owner());
      }
    }
  }
}

//-----------------------------------------------------------------------------
void UniformRefinement::operator()(Mesh& mesh)
{
  // Create new mesh, refinement manager and open for editing
  Mesh refined_mesh(mesh.type(), mesh.space());
  MeshEditor editor(refined_mesh, refined_mesh.type(), refined_mesh.space());

  // Refinement pattern provides the number of refined vertices
  editor.init_vertices(mesh.type().RefinementPattern::num_refined_vertices(mesh));

  // Refinement pattern provides the number of refined cells
  editor.init_cells(mesh.type().RefinementPattern::num_refined_cells(mesh));

  // Add vertices for each entity, skip cell-based vertices for point clouds
  add_refined_vertices<Vertex>(editor, mesh);
  add_refined_vertices<Edge>  (editor, mesh);
  add_refined_vertices<Face>  (editor, mesh);
  if (mesh.type().dim() > 0) { add_refined_vertices<Cell>  (editor, mesh); }

  // Add cells
  uint current_cell = 0;
  for (Cell::iterator c(mesh); !c.end(); ++c)
  {
    mesh.type().refine_cell(*c, editor, current_cell);
  }

  // Apply numbering of new entities and close edition
  editor.close();

  // Re-map local vertices contiguously
  uint const num_local = refined_mesh.topology().size(0);
  uint vertex_count = 0;
  Array<uint> vertex_map(num_local, num_local);
  for (Cell::iterator c(refined_mesh); !c.end(); ++c)
  {
    for (Vertex::iterator v(*c); !v.end(); ++v)
    {
      dolfin_assert(v->index() < num_local);
      if (vertex_map[v->index()] == num_local)
      {
        vertex_map[v->index()] = vertex_count++;
      }
    }
  }

  // Reorder geometry
  dolfin_assert(vertex_count == refined_mesh.geometry().size());
  refined_mesh.geometry().remap(vertex_map);

  // Reorder connectivities
  dolfin_assert(vertex_count == refined_mesh.topology().size(0));
  refined_mesh.topology().remap(0, vertex_map);

  // Overwrite old mesh with refined mesh
  swap( mesh, refined_mesh );
  mesh.topology().renumber();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

