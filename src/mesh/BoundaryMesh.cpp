// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <iostream>

#include <dolfin/main/PE.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
BoundaryMesh::BoundaryMesh(Mesh& mesh, BoundaryMesh::Type type) :
    Mesh(mesh.type(), mesh.space(), mesh.topology().comm()),
    MeshDependent(static_cast<Mesh&>(mesh)),
    type_(type),
    boundary_of_boundary_(false),
    cell_map_(),
    vertex_map_(),
    subdomain_(NULL)
{
  init(mesh, type);
}
//-----------------------------------------------------------------------------
BoundaryMesh::BoundaryMesh(BoundaryMesh& mesh, BoundaryMesh::Type type) :
    Mesh(mesh.type(), mesh.space(), mesh.topology().comm()),
    MeshDependent(static_cast<Mesh&>(mesh)),
    type_(type),
    boundary_of_boundary_(true),
    cell_map_(),
    vertex_map_(),
    subdomain_(NULL)
{
  init(mesh, type);
}
//-----------------------------------------------------------------------------
BoundaryMesh::BoundaryMesh(Mesh& mesh, SubDomain const& subdomain,
                           BoundaryMesh::Type type) :
    Mesh(mesh.type(), mesh.space(), mesh.topology().comm()),
    MeshDependent(mesh),
    type_(type),
    boundary_of_boundary_(false),
    cell_map_(),
    vertex_map_(),
    subdomain_(&subdomain)
{
  init(mesh, type);
}
//-----------------------------------------------------------------------------
BoundaryMesh::BoundaryMesh(BoundaryMesh const& other) :
    Mesh(other),
    MeshDependent(other.mesh()),
    type_(other.type_),
    boundary_of_boundary_(other.boundary_of_boundary_),
    cell_map_(other.cell_map_),
    vertex_map_(other.vertex_map_),
    subdomain_(other.subdomain_)
{
}
//-----------------------------------------------------------------------------
void BoundaryMesh::init(Mesh& mesh, BoundaryMesh::Type type)
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
BoundaryMesh::BoundaryMesh(BoundaryMesh& boundary, SubDomain const& subdomain,
                           bool inside) :
    Mesh(boundary.type(), boundary.space(), boundary.topology().comm()),
    MeshDependent(boundary.mesh()),
    type_(boundary.type_),
    boundary_of_boundary_(true),
    cell_map_(),
    vertex_map_(),
    subdomain_(&subdomain)
{
  Mesh& mesh = boundary.mesh();
  uint const gdim = mesh.geometry_dimension();
  uint const tdim = mesh.topology_dimension();

  if (tdim == 1)
  {
    vertex_map_.clear();
    for (VertexIterator v(boundary); !v.end(); ++v)
    {
      if (subdomain_->inside(v->x(), !v->is_shared()))
      {
        vertex_map_.push_back(boundary.vertex_map_[v->index()]);
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
    if (mesh.is_distributed())
    {
      this->distdata()[0].assign(mesh.distdata()[0], vertex_map_);
    }
    editor.close();
  }
  else
  {

    uint const num_verts = boundary.vertex_map_.size();
    Array<uint> boundary_vertices(num_verts, num_verts);
    for (CellIterator c(boundary); !c.end(); ++c)
    {
      bool const on_boundary = !c->is_shared();
      if ((inside && subdomain_->enclosed(*c, on_boundary))
          || (!inside && subdomain_->overlap(*c, on_boundary)))
      {
        cell_map_.push_back(c->index());
        for (VertexIterator v(*c); !v.end(); ++v)
        {
          uint const vertex_index = v->index();
          if (boundary_vertices[v->index()] == num_verts)
          {
            boundary_vertices[vertex_index] = vertex_map_.size();
            vertex_map_.push_back(boundary.vertex_map_[vertex_index]);
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
    editor.init_cells(cell_map_.size());
    //
    uint const d = mesh.type().facet_dim();
    uint const num_facet_vertices = mesh.type().num_vertices(d);
    uint * facet_vertices = new uint[num_facet_vertices];
    for (uint i = 0; i < cell_map_.size(); ++i)
    {
      Cell c(boundary, cell_map_[i]);
      cell_map_[i] = boundary.cell_map_[cell_map_[i]];
      for (VertexIterator v(c); !v.end(); ++v)
      {
        facet_vertices[v.pos()] = boundary_vertices[v->index()];
      }
      editor.add_cell(i, &facet_vertices[0]);
    }
    delete[] facet_vertices;

    // If the mesh is distributed, set global numbering and copy the ownership
    if (mesh.is_distributed())
    {
      this->distdata()[0].assign(mesh.distdata()[0], vertex_map_);
      this->distdata()[d].assign(mesh.distdata()[d], cell_map_);
    }
    editor.close();
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
bool BoundaryMesh::is_boundary_of_boundary() const
{
  return boundary_of_boundary_;
}
//-----------------------------------------------------------------------------
void BoundaryMesh::compute(Mesh& mesh, bool exterior, bool interior)
{
  uint const gdim = mesh.geometry_dimension();
  uint const tdim = mesh.topology_dimension();

  // If the boundary is full then no need to compute the facet map
  bool const full = interior && exterior;

  message(1, "BoundaryMesh : compute %s boundary",
          (full ? "full" : (interior ? "interior" : "exterior")));

  if (tdim <= 1)
  {
    vertex_map_.clear();
    for (VertexIterator v(mesh); !v.end(); ++v)
    {
      // Boundary facets are connected to exactly one cell
      bool const bndr = (v->num_entities(tdim) == 1);
      // Interior facets are shared while exterior facets are not
      bool const shrd = v->is_shared();
      if (bndr && (full || (interior && shrd) || (exterior && !shrd)))
      {
        if ((subdomain_ == NULL) || (subdomain_->inside(v->x(), bndr && !shrd)))
        {
          vertex_map_.push_back(v->index());
        }
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

    uint const pe_size = PE::size();
    Array<uint> * shared_vertices = new Array<uint>[pe_size];
    Array<uint> boundary_vertices(num_verts, num_verts);
    for (FacetIterator f(mesh); !f.end(); ++f)
    {
      // Boundary facets are connected to exactly one cell
      bool const bndr = (f->num_entities(tdim) == 1);
      // Interior facets are shared while exterior facets are not
      bool const shrd = f->is_shared();
#if DEBUG
      if(shrd && !f->has_all_vertices_shared())
      {
        error("Shared facets have not all vertices shared");
      }
#endif
      if (bndr && (full || (interior && shrd) || (exterior && !shrd)))
      {
        if (subdomain_ != NULL)
        {
          bool inside = false;
          for (VertexIterator v(*f); !v.end(); ++v)
          {
            if(subdomain_->inside(v->x(), bndr && !shrd))
            {
              inside = true;
            }
          }
          if(!inside)
          {
            continue;
          }
        }
        cell_map_.push_back(f->index());
        for (VertexIterator v(*f); !v.end(); ++v)
        {
          uint const vertex_index = v->index();
          if (boundary_vertices[vertex_index] == num_verts)
          {
            boundary_vertices[vertex_index] = vertex_map_.size();
            vertex_map_.push_back(vertex_index);
            // Collect shared vertices per adjacent ranks
            if (v->is_shared())
            {
              _set<uint> const& adjs = mesh.distdata()[0].get_shared_adj(vertex_index);
              for (_set<uint>::const_iterator adj = adjs.begin(); adj != adjs.end();
                   ++adj)
              {
                dolfin_assert(mesh.distdata()[0].get_adj_ranks().count(*adj));
                shared_vertices[*adj].push_back(v->global_index());
              }
            }
          }
        }
      }
    }

    // The mesh distribution is not constrained: non-owner of exterior facets
    // may own lower-dimensional entities on the exterior boundary.
    if(mesh.is_distributed())
    {
#if HAVE_MPI
      MPI_Status status;
      DistributedData const& distdata = mesh.distdata()[0];
      _set<uint> const& vadjs = distdata.get_adj_ranks();
      uint recvmax = shared_vertices[0].size();
      for (uint j = 1; j < pe_size; ++j)
      {
        recvmax = std::max(recvmax, (uint) shared_vertices[j].size());
      }
      MPI::all_reduce<MPI::max>(recvmax, recvmax);
      uint * recvbuf = new uint[recvmax];
      int recvcount;

      _set<uint> added_vertices(vertex_map_.begin(), vertex_map_.end());

      // If one adjacent rank has no boundary cell but the current rank has
      // ghost entities owned by this rank then any algorithm based on mesh
      // facets is bound to fail.
      // Naive implementation for now, just send all the ghost vertices
      for (_set<uint>::const_iterator adj = vadjs.begin(); adj != vadjs.end();
           ++adj)
      {
        MPI::check_error( MPI_Send(&shared_vertices[(*adj)][0],
                                   shared_vertices[(*adj)].size(),
                                   MPI_UNSIGNED, (*adj), 0, MPI::DOLFIN_COMM) );
      }
      for (_set<uint>::const_iterator adj = vadjs.begin(); adj != vadjs.end();
           ++adj)
      {
        MPI::check_error( MPI_Recv(&recvbuf[0], recvmax, MPI_UNSIGNED, (*adj),
                                   0, MPI::DOLFIN_COMM, &status) );
        MPI::check_error( MPI_Get_count(&status, MPI_UNSIGNED, &recvcount) );
        for(int k = 0; k < recvcount; ++k)
        {
          dolfin_assert(distdata.has_global(recvbuf[k]));
          uint const local_index = distdata.get_local(recvbuf[k]);
          if(added_vertices.count(local_index) == 0)
          {
            vertex_map_.push_back(local_index);
            added_vertices.insert(local_index);
          }
        }
      }
      //
      delete [] recvbuf;

#endif /* HAVE_MPI */
    }

    delete [] shared_vertices;

    // Create boundary vertices and cells
    MeshEditor editor(*this, mesh.type().facetType(), gdim);
    editor.init_vertices(vertex_map_.size());
    for (uint i = 0; i < vertex_map_.size(); ++i)
    {
      editor.add_vertex(i, mesh.geometry().x(vertex_map_[i]));
    }
    editor.init_cells(cell_map_.size());
    uint const d = mesh.type().facet_dim();
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
  }

  message(1, "BoundaryMesh : number of cells = %u, number of vertices %u",
          cell_map_.size(), vertex_map_.size());
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
