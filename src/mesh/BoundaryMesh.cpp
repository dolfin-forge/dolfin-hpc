// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson 2008.
//
// First added:  2006-06-21
// Last changed: 2008-06-26

#include <iostream>

#include <dolfin/log/log.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
BoundaryMesh::BoundaryMesh(Mesh& mesh, BoundaryMesh::Type type) :
    Mesh(),
    MeshDependent(mesh),
    type_(type),
    cell_map_(),
    vertex_map_()
{
  switch (type)
    {
    case BoundaryMesh::exterior:
      // Exterior boundary i.e facets at the domain boundary
      compute(mesh, true, false);
      break;
    case BoundaryMesh::interior:
      // Interior boundary i.e facets between processors
      compute(mesh, false, true);
      break;
    case BoundaryMesh::full:
      // Full boundary including facets between processors
      compute(mesh, true, true);
      break;
    default:
      error("Unknown boundary mesh type.");
      break;
    }
}
//-----------------------------------------------------------------------------
BoundaryMesh::~BoundaryMesh()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
uint BoundaryMesh::facet_index(Cell const& boundary_cell) const
{
  dolfin_assert(&boundary_cell.mesh() == this);
  return cell_map_[boundary_cell.index()];
}
//-----------------------------------------------------------------------------
uint BoundaryMesh::facet_index(uint boundary_cell_index) const
{
  dolfin_assert(boundary_cell_index < cell_map_.size());
  return cell_map_[boundary_cell_index];
}
//-----------------------------------------------------------------------------
uint BoundaryMesh::vertex_index(Vertex const& boundary_vertex) const
{
  dolfin_assert(&boundary_vertex.mesh() == this);
  return vertex_map_[boundary_vertex.index()];
}
//-----------------------------------------------------------------------------
uint BoundaryMesh::vertex_index(uint boundary_vertex_index) const
{
  dolfin_assert(boundary_vertex_index < vertex_map_.size());
  return vertex_map_[boundary_vertex_index];
}
//-----------------------------------------------------------------------------
BoundaryMesh::Type BoundaryMesh::boundary_type() const
{
  return type_;
}
//-----------------------------------------------------------------------------
void BoundaryMesh::compute(Mesh& mesh, bool exterior, bool interior)
{
  uint const gdim = mesh.geometry().dim();
  uint const tdim = mesh.topology().dim();

  // If the boundary is full then no need to compute the facet map
  bool const full = interior && exterior;

  message(1, "BoundaryMesh : compute %s boundary",
          (full ? "full" : (interior ? "interior" : "exterior")));

  if (tdim == 1)
  {
    vertex_map_.clear();
    for (VertexIterator v(mesh); !v.end(); ++v)
    {
      // Boundary facets are connected to exactly one cell
      if ((full && (v->num_entities(tdim) == 1)) ||
          (!full&& ((interior && v->is_shared()) ||
                    (exterior && !v->is_shared()))))
      {
        vertex_map_.push_back(v->index());
      }
    }
    cell_map_ = vertex_map_;

    // Create boundary vertices and cells
    MeshEditor editor(*this, mesh.type().facetType(), gdim);
    editor.init_vertices(vertex_map_.size());
    editor.init_cells(cell_map_.size());
    MeshGeometry const& geom = mesh.geometry();
    for (uint i = 0; i < vertex_map_.size(); ++i)
    {
      editor.add_vertex(i, geom.x(vertex_map_[i]));
      editor.add_cell(i, &vertex_map_[i]);
    }
    // If the mesh is distributed, set global numbering and copy the ownership
    if(mesh.is_distributed())
    {
      this->distdata()[0].assign(mesh.distdata()[0], vertex_map_);
    }
    editor.close();
  }
  else
  {
    uint const num_verts = mesh.size(0);
    cell_map_.clear();
    vertex_map_.clear();

    Array<uint> boundary_vertices(num_verts, num_verts);
    for (FacetIterator f(mesh); !f.end(); ++f)
    {
      // Boundary facets are connected to exactly one cell
      if ((f->num_entities(tdim) == 1) &&
          (full || ((interior && f->is_shared()) || (exterior && !f->is_shared()))))
      {
        cell_map_.push_back(f->index());
        for (VertexIterator v(*f); !v.end(); ++v)
        {
          uint const vertex_index = v->index();
          if (boundary_vertices[vertex_index] == num_verts)
          {
            boundary_vertices[vertex_index] = vertex_map_.size();
            vertex_map_.push_back(vertex_index);
          }
        }
      }
    }

    // Create boundary vertices and cells
    MeshEditor editor(*this, mesh.type().facetType(), gdim);
    editor.init_vertices(vertex_map_.size());
    for (uint i = 0; i < vertex_map_.size(); ++i)
    {
      editor.add_vertex(i, mesh.geometry().x(vertex_map_[i]));
    }
    editor.init_cells(cell_map_.size());;
    uint const d = tdim - 1;
    uint const num_facet_vertices = mesh.type().num_vertices(d);
    uint * facet_vertices = new uint[num_facet_vertices];
    CellType const& celltype = mesh.type();
    for (uint i = 0; i < cell_map_.size(); ++i)
    {
      Facet facet(mesh, cell_map_[i]);
      for (uint v = 0; v < num_facet_vertices; ++v)
      {
        facet_vertices[v] = boundary_vertices[facet.entities(0)[v]];
      }
      // Reorder vertices so facet is right-oriented w.r.t. facet normal
      celltype.order_facet(&facet_vertices[0], facet);
      editor.add_cell(i, &facet_vertices[0]);
    }
    delete [] facet_vertices;
    // If the mesh is distributed, set global numbering and copy the ownership
    if(mesh.is_distributed())
    {
      this->distdata()[0].assign(mesh.distdata()[0], vertex_map_);
      this->distdata()[d].assign(mesh.distdata()[d], cell_map_);
    }
    editor.close();
    boundary_vertices.clear();
  }

  if(mesh.is_distributed())
  {
#if HAVE_MPI
    _set<uint> adjs = mesh.distdata()[0].get_adj_ranks();
    MPI_Status status;
    uint sendbuf = cell_map_.size();
    for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it)
    {
      uint recvbuf;
      MPI_Sendrecv(&sendbuf, 1, MPI_UNSIGNED, *it, 0,
                   &recvbuf, 1, MPI_UNSIGNED, *it, 0,
                   MPI::DOLFIN_COMM, &status);
      if (recvbuf == 0)
      {
        error("BoundaryMesh : adjacent rank %u has no cell, case unimplemented",
              *it);
      }
    }
#endif /* HAVE_MPI */

  }

  message(1, "BoundaryMesh : number of cells = %u, number of vertices %u",
          cell_map_.size(), vertex_map_.size());
}
//-----------------------------------------------------------------------------
void BoundaryMesh::check() const
{
  Mesh& mesh = const_cast<Mesh&>(this->mesh());
  uint const mdim = mesh.topology().dim();
  uint const fdim = (mdim > 0 ? mdim - 1 : 0);
  uint const bdim = this->topology().dim();
  if (bdim != fdim)
  {
    error("BoundaryMesh : invalid dimension %u for boundary mesh, expected %u.",
          bdim, fdim);
  }

  // Check interior entities
  if(type_ == BoundaryMesh::interior || type_ == BoundaryMesh::full)
  {
    for (uint d = 0; d <= bdim; ++d)
    {
      this->check_interior(d);
    }
  }

  // Check exterior entities
  if (type_ == BoundaryMesh::exterior || type_ == BoundaryMesh::full)
  {
    for (uint d = 0; d <= bdim; ++d)
    {
      this->check_exterior(d);
    }
  }

}
//-----------------------------------------------------------------------------
void BoundaryMesh::check_interior(uint dim) const
{
  message("BoundaryMesh : check interior boundary entities of dimension %d",
          dim);

  /**
   *  CHECK:
   *
   *  The interior boundary consists of the facets shared between processes and
   *  cannot be empty for a parallel run (provided that the domain covered by
   *  the mesh is made of one piece).
   *
   */

  Mesh& mesh = const_cast<Mesh&>(this->mesh());
  uint const mdim = mesh.topology().dim();
  uint const bdim = this->topology().dim();

  //
  if (dim > bdim)
  {
    error("BoundaryMesh : invalid dimension %u for interior check.", dim);
  }
  //
  if (this->size(bdim) == 0 && mesh.is_distributed())
  {
    error("BoundaryMesh : distributed mesh has no interior boundary.");
  }
  //
  if (this->size(bdim) > 0 && !mesh.is_distributed())
  {
    error("BoundaryMesh : non-distributed mesh has an interior boundary.");
  }
  //
  if(this->topology().num_shared(dim) != this->topology().size(dim))
  {
    error("Inconsistent number of entities: (shared) %d  != %d (boundary)",
          this->topology().num_shared(dim), this->topology().size(dim));
  }
  // Test all the mesh entities at the interior boundary
  // All the entities should be shared and some are ghosted
  if(mesh.is_distributed())
  {
    uint num_invalid_facets = 0;
    uint num_invalid_shared = 0;
    if (dim == bdim)
    {
      for (uint c = 0; c < this->size(bdim); ++c)
      {
        Facet facet(mesh, this->facet_index(c));
        if (!facet.is_shared())
        {
          ++num_invalid_shared;
        }
        if (facet.num_entities(mdim) != 1)
        {
          ++num_invalid_facets;
        }
      }
    }
    else
    {
      for (uint c = 0; c < this->size(bdim); ++c)
      {
        Facet facet(mesh, this->facet_index(c));
        for (MeshEntityIterator e(facet, dim); !e.end(); ++e)
        {
          if (!e->is_shared())
          {
            ++num_invalid_shared;
          }
        }
      }
    }
    //
    if (num_invalid_shared > 0)
    {
      error("Interior boundary entities of dim %d: %d invalid shared data.",
            dim, num_invalid_shared);
    }
    if (num_invalid_facets > 0)
    {
      error("Interior boundary entities of dim %d: %d invalid connectivity.",
            dim, num_invalid_facets);
    }
  }
}
//-----------------------------------------------------------------------------
void BoundaryMesh::check_exterior(uint dim) const
{
  message("BoundaryMesh : check exterior boundary entities of dimension %d",
          dim);

  /**
   *  CHECK:
   *
   *  The exterior boundary consists of the facets located on the boundary of
   *  the domain and thus not shared between processes.
   *
   */

  Mesh& mesh = const_cast<Mesh&>(this->mesh());
  uint const mdim = mesh.topology().dim();
  uint const bdim = this->topology().dim();

  uint num_invalid_facets = 0;
  uint num_invalid_shared = 0;
  if (dim == bdim)
  {
    for (uint c = 0; c < this->size(bdim); ++c)
    {
      Facet facet(mesh, this->facet_index(c));
      if (facet.is_shared())
      {
        ++num_invalid_shared;
      }
      if (facet.num_entities(mdim) != 1)
      {
        ++num_invalid_facets;
      }
    }
  }
  else
  {
    warning("BoundaryMesh : unimplemented for dimension %u", dim);
  }
  //
  if (num_invalid_shared > 0)
  {
    error("Interior boundary entities of dim %d: %d invalid shared data.",
          dim, num_invalid_shared);
  }
  if (num_invalid_facets > 0)
  {
    error("Interior boundary entities of dim %d: %d invalid connectivity.",
          dim, num_invalid_facets);
  }
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
