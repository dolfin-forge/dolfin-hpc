// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
// Modified by Stefanie Strunk, 2013.
// Modified by Aurelien Larcher, 2017.
//
// First added:  2006-06-08
// Last changed: 2017-12-15

#include <dolfin/mesh/UniformMeshRefinement.h>

#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Face.h>
#include <dolfin/mesh/Cell.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
template<class E>
void add_refined_vertices(MeshEditor& editor, Mesh& mesh)
{
  Mesh& refined_mesh = editor.mesh();
  uint const tdim = refined_mesh.topology().dim();
  uint const edim = entity_dimension<E>(mesh);
  if (tdim > edim && mesh.type().refinement_needs_entities(edim))
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
      for (typename E::shared it(mesh); !it.end(); ++it)
      {
        dist.setall_shared_adj(voffset + it.index(), it.adj());
      }
      for (typename E::ghost it(mesh); !it.end(); ++it)
      {
        dist.set_ghost(voffset + it.index(), it.owner());
      }
    }
  }
}

//-----------------------------------------------------------------------------
void UniformMeshRefinement::refine(Mesh& mesh)
{
  message(1, "Refining %s mesh uniformly.", mesh.type().str().c_str());

  // Create new mesh, refinement manager and open for editing
  Mesh refined_mesh;
  MeshEditor editor(refined_mesh, mesh.type(), mesh.space());

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

  // Overwrite old mesh with refined mesh
  mesh.swap(refined_mesh);
  mesh.topology().renumber();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

