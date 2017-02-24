// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-10-05
// Last changed: 2014-10-05

#include <dolfin/mesh/MappedManifold.h>

#include <dolfin/mesh/PeriodicSubDomain.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MappedManifold::MappedManifold(Mesh& mesh, PeriodicSubDomain const& subdomain) :
    Mesh(),
    MeshDependent(mesh),
    subdomainG_(subdomain)
{
  init();
}

//-----------------------------------------------------------------------------
MappedManifold::~MappedManifold()
{
}

//-----------------------------------------------------------------------------
uint MappedManifold::facet_index(Cell const& boundary_cell)
{
  dolfin_assert(&boundary_cell.mesh() == this);
  return cell_map_[boundary_cell.index()];
}

//-----------------------------------------------------------------------------
uint MappedManifold::vertex_index(Vertex const& boundary_vertex)
{
  dolfin_assert(&boundary_vertex.mesh() == this);
  return vertex_map_[boundary_vertex.index()];
}

//-----------------------------------------------------------------------------
void MappedManifold::init()
{
  //
  Mesh& globalmesh = this->mesh();
  BoundaryMesh& boundary = globalmesh.exterior_boundary();

  uint const tdim = globalmesh.topology().dim();
  uint const gdim = globalmesh.geometry().dim();
  MeshEditor editor(*this, boundary.type().cellType(),
                    boundary.geometry().dim());
  Array<uint> mm_vertices(boundary.size(0));
  uint const invalid_vertex_index = mm_vertices.size();
  mm_vertices = invalid_vertex_index;
  Array<real> mm_coordinates(gdim * mm_vertices.size());
  uint const num_cell_vertices = boundary.type().num_entities(0);
  Array<uint> mm_cell_vertices(boundary.num_cells() * num_cell_vertices);

  Array<uint> localGH;

  //
  uint mm_cell_count = 0;
  uint mm_vertex_count = 0;
  Point x1;
  cell_map_.clear();
  for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
  {
    Facet facet(globalmesh, boundary.facet_index(*bcell));
    Cell cell(globalmesh, facet.entities(tdim)[0]);

    //
    bool is_Gfacet = false;
    bool is_Hfacet = false;
    uint Hfound = 0;
    for (VertexIterator v(*bcell); !v.end(); ++v)
    {
      is_Gfacet |= subdomainG_.inside(v->x(), true);
      subdomainG_.map(v->x(), &x1[0]);
      if(subdomainG_.inside(&x1[0], true))
      {
        ++Hfound;
      }
    }
    //is_Hfacet = (Hfound == facet.numEntities(0));
    is_Hfacet = (Hfound > 0);

    if (is_Gfacet)
    {
      facetsG_.insert(facet.index());
    }

    if (is_Hfacet)
    {
      facetsH_.insert(facet.index());

      // Add cell (boundary facet)
      cell_map_.push_back(facet.index());
      // Add vertices expressed in the reference manifold coordinates
      for (VertexIterator v(*bcell); !v.end(); ++v)
      {
        if (mm_vertices[v->index()] == invalid_vertex_index)
        {
          mm_vertices[v->index()] = mm_vertex_count;
          subdomainG_.map(v->x(), &mm_coordinates[mm_vertex_count * gdim]);
          ++mm_vertex_count;

          Array<uint> matching_facets;
          boundary.intersector().overlap(x1, matching_facets);
          if(!matching_facets.empty())
          {
            // Local H -> G mapping
            facetsL_.insert(facet.index());
          }
        }
        mm_cell_vertices[mm_cell_count * num_cell_vertices + v.pos()] =
            mm_vertices[v->index()];
      }
      ++mm_cell_count;
    }

    if (is_Gfacet && is_Hfacet)
    {
      facetsI_.insert(facet.index());
    }
  }

  //
  editor.init_vertices(mm_vertex_count);
  dolfin_assert(cell_map_.size() == mm_cell_count);
  editor.init_cells(mm_cell_count);

  // Initialize mapping from vertices in boundary to vertices in mesh
  vertex_map_.resize(mm_vertex_count);

  // Create vertices
  for (uint v = 0; v < mm_vertices.size(); ++v)
  {
    uint const vertex_index = mm_vertices[v];
    if (vertex_index != invalid_vertex_index)
    {
      // Create mapping from boundary vertex to mesh vertex
      vertex_map_[vertex_index] = v;

      // Add vertex with coordinates expressed in reference manifold
      editor.add_vertex(vertex_index, &mm_coordinates[vertex_index * gdim]);
    }
  }

  // Create cells (boundary mesh is already ordered)
  Array<uint> cell_vertices(num_cell_vertices);
  for (uint c = 0; c < mm_cell_count; ++c)
  {
    // Add cell
    editor.add_cell(c, &mm_cell_vertices[c * num_cell_vertices]);
  }
  editor.close();

  //
  message(1, "Computed reference manifold with %d cell and %d vertices.",
          mm_cell_count, mm_vertex_count);
  message(1, "Generated sets of %8d G facets, %8d H facets, %8d I facet",
          facetsG_.size(), facetsH_.size(), facetsI_.size());
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
