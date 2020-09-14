// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/QuadrilateralCell.h>

#include <dolfin/common/maybe_unused.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>

#include <algorithm>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
uint const QuadrilateralCell::NE[3][3] =
{ { 1, 0, 0 }, { 2, 1, 0 }, { 4, 4, 1 } };

// UFC: Vertex Coordinates
real const QuadrilateralCell::VC[4][2] =
{ { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };

// UFC: Edge - Incident Vertices
uint const QuadrilateralCell::EIV[4][2] =
{ { 2, 3 }, { 1, 2 }, { 0, 3 }, { 0, 1 } };

// UFC: Edge - Non-Incident Vertices
uint const QuadrilateralCell::ENV[4][2] =
{ { 0, 1 }, { 0, 3 }, { 1, 2 }, { 2, 3 } };

//-----------------------------------------------------------------------------
QuadrilateralCell::QuadrilateralCell() :
    CellType("quadrilateral", CellType::quadrilateral, CellType::interval)
{
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::create_entities(uint** e, uint dim, uint const* v) const
{
  // We only need to know how to create edges
  if (dim != 1)
  {
    error("Invalid topological dimension for creation of entities: %d.", dim);
  }

  // Create the four edges using UFC ordering
  e[0][0] = v[2];
  e[0][1] = v[3];
  e[1][0] = v[1];
  e[1][1] = v[2];
  e[2][0] = v[0];
  e[2][1] = v[3];
  e[3][0] = v[0];
  e[3][1] = v[1];
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::order_entities(MeshTopology& topology, uint i) const
{
  // Sort i - j for i > j: 1 - 0, 2 - 0, 2 - 1
  dolfin_assert(topology.type().cellType() == this->cell_type);

  // Sort local vertices on edges in ascending order, connectivity 1 - 0
  if (topology.connectivity(1, 0))
  {
    dolfin_assert(topology.connectivity(2, 1));

    // Get edges
    Array<uint> const & cell_edges = topology(2, 1)[i];

    // Sort vertices on each edge
    for (uint i = 0; i < 4; ++i)
    {
      Array<uint> & edge_vertices = topology(1, 0)[cell_edges[i]];
      std::sort(edge_vertices.data(), edge_vertices.data() + 2);
    }
  }

  // Vertices are supposed to have be ordered appropriately to obtain a
  // conforming quadrilateral

  // Sort local edges on cell after non-incident vertex, connectivity 2 - 1
  if (topology.connectivity(2, 1))
  {
    dolfin_assert(topology.connectivity(2, 1));

    // Get cell vertices and edges
    Array<uint> const & cell_vertices = topology(2, 0)[i];
    Array<uint> & cell_edges = topology(2, 1)[i];

    // Loop on non-incident vertices on cell as lexicographical pairs
    // (i, j): (0,1) (0,3) (1,2) (2,3)
    uint m = 0;
    for (uint t = 0; t < 4; ++t)
    {
      uint v0 = ENV[t][0];
      uint v1 = ENV[t][1];

      // Loop on edges
      for (uint k = m; k < 4; ++k)
      {
        // Get local vertices on edge
        Array<uint> const & edge_vertices = topology(1, 0)[cell_edges[k]];

        // Check if the ith and jth vertex of the cell are non-incident on edge k
        if (!std::count(edge_vertices.data(), edge_vertices.data() + 2,
                        cell_vertices[v0])
         && !std::count(edge_vertices.data(), edge_vertices.data() + 2,
                        cell_vertices[v1]))
        {
          // Swap edge numbers
          uint tmp = cell_edges[m];
          cell_edges[m] = cell_edges[k];
          cell_edges[k] = tmp;
          ++m;
          break;
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::order_facet(uint vertices[], Facet& facet) const
{
  // Get mesh
  Mesh& mesh = facet.mesh();

  // Get the vertex opposite to the facet (the one we remove)
  uint vertex = 0;
  const Cell cell(mesh, facet.entities(mesh.topology_dimension())[0]);
  for (uint i = 0; i < cell.num_entities(0); i++)
  {
    bool not_in_facet = true;
    vertex = cell.entities(0)[i];
    for (uint j = 0; j < facet.num_entities(0); j++)
    {
      if (vertex == facet.entities(0)[j])
      {
        not_in_facet = false;
        break;
      }
    }
    if (not_in_facet)
    {
      break;
    }
  }

  // Order
  Point const p = mesh.geometry().point(vertex);
  Point p0 = mesh.geometry().point(facet.entities(0)[0]);
  Point p1 = mesh.geometry().point(facet.entities(0)[1]);
  Point v = p1 - p0;
  Point n(v[1], -v[0]);
  if (n.dot(p0 - p) < 0.0)
  {
    std::swap(vertices[0], vertices[1]);
  }
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::refine_cell(Cell& cell, MeshEditor& editor,
                                   uint& current_cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get vertices and edges
  Array<uint> const & v = cell.entities(0);
  dolfin_assert(!v.empty());
  Array<uint> const & e = cell.entities(1);
  dolfin_assert(!e.empty());

  // Compute indices for the nine new vertices
  uint const v0 = v[0];
  uint const v1 = v[1];
  uint const v2 = v[2];
  uint const v3 = v[3];
  uint const eoffset = cell.mesh().size(0);
  uint const e0 = eoffset + e[findEdge(0, cell)];
  uint const e1 = eoffset + e[findEdge(1, cell)];
  uint const e2 = eoffset + e[findEdge(2, cell)];
  uint const e3 = eoffset + e[findEdge(3, cell)];
  uint const coffset = eoffset + cell.mesh().size(1);
  uint const c0 = coffset + cell.index();

  // Add the four new cells
  uint cv0[4] = { v0, e3, c0, e2 };
  editor.add_cell(current_cell++, &cv0[0]);
  uint cv1[4] = { e3, v1, e1, c0};
  editor.add_cell(current_cell++, &cv1[0]);
  uint cv2[4] = { c0, e1, v2, e0 };
  editor.add_cell(current_cell++, &cv2[0]);
  uint cv3[4] = { e2, c0, e0, v3 };
  editor.add_cell(current_cell++, &cv3[0]);
}
//-----------------------------------------------------------------------------
bool QuadrilateralCell::intersects(MeshEntity const& e, Point const&) const
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.num_entities(0) == NE[2][0]);
  MAYBE_UNUSED(e);

  // Get the coordinates of the vertices
  /*
  MeshGeometry const& geometry = e.mesh().geometry();
  uint const* vertices = e.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
   */

  error("Collision of quadrilateral with point no implemented.");

  return true;
}
//-----------------------------------------------------------------------------
bool QuadrilateralCell::intersects(MeshEntity const& e, Point const&,
                                   Point const&) const
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.num_entities(0) == NE[2][0]);
  MAYBE_UNUSED(e);

  // Get the coordinates of the vertices
  /*
  MeshGeometry const& geometry = e.mesh().geometry();
  uint const* vertices = e.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
   */

  error("Collision of quadrilateral with segment not implemented");

  return true;
}
//-----------------------------------------------------------------------------
std::string QuadrilateralCell::description() const
{
  return std::string("quadrilateral (hypercube of topological dimension 2)");
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::create_reference_cell(Mesh& mesh) const
{
  MeshEditor me(mesh, CellType::quadrilateral, 2, DOLFIN_COMM_SELF);
  me.init_vertices(4);
  me.add_vertex(0, VC[0]);
  me.add_vertex(1, VC[1]);
  me.add_vertex(2, VC[2]);
  me.add_vertex(3, VC[3]);
  me.init_cells(1);
  uint const cv0[4] = { 0, 1, 2, 3 };
  me.add_cell(0, cv0);
  me.close();
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::disp() const
{
  message("QuadrilateralCell");
  begin(  "-----------------");
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
bool QuadrilateralCell::check(Cell& cell) const
{
  bool ret = CellType::check(cell);

  // Check edge -> incident vertices mapping
  if (cell.mesh().topology().connectivity(1, 0))
  {
    Array<uint> const & v = cell.entities(0);
    dolfin_assert( not v.empty() );
    Array<uint> const & e = cell.entities(1);
    dolfin_assert( not e.empty() );
    MeshTopology const& topology = cell.mesh().topology();
    for (uint i = 0; i < 4; ++i)
    {
      Array<uint> const & ev = topology(1, 0)[e[i]];
      dolfin_assert(ev.size() >= 2);
      for (uint j = 0; j < 2; ++j)
      {
        if (ev[j] != v[EIV[i][0]] && ev[j] != v[EIV[i][1]])
        {
          ret = false;
          warning("CellType : invalid edge -> incident vertices mapping\n"
                "e[%u] = (v%u, v%u)", i, ev[0], ev[1]);
        }
      }
    }
  }

  return ret;
}
//-----------------------------------------------------------------------------
uint QuadrilateralCell::findEdge(uint i, Cell const& cell) const
{
  // Ordering convention for edges (order of non-incident vertices)

  // Get vertices and edges
  Array<uint> const & v = cell.entities(0);
  dolfin_assert( not v.empty() );
  Array<uint> const & e = cell.entities(1);
  dolfin_assert( not e.empty() );

  // Look for edge satisfying ordering convention
  MeshTopology const& topology = cell.mesh().topology();
  uint const v0 = v[ENV[i][0]];
  uint const v1 = v[ENV[i][1]];
  for (uint j = 0; j < 4; ++j)
  {
    Array<uint> const & ev = topology(1, 0)[e[j]];
    dolfin_assert( not ev.empty() );
    if (ev[0] != v0 && ev[0] != v1 && ev[1] != v0 && ev[1] != v1)
    {
      return j;
    }
  }

  // We should not reach this
  error("Unable to find edge with index %d in quadrilateral.", cell.index());

  return 0;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
