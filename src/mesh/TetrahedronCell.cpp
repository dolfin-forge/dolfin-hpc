// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/TetrahedronCell.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/maybe_unused.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/GeometricPredicates.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>

#include <algorithm>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
uint const TetrahedronCell::NE[4][4] =
{ { 1, 0, 0, 0 }, { 2, 1, 0, 0 }, { 3, 3, 1, 0 }, { 4, 6, 4, 1 } };

// UFC: Vertex Coordinates
real const TetrahedronCell::VC[4][3] =
{ { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

// UFC: Edge - Incident Vertices
uint const TetrahedronCell::EIV[6][2] =
{ { 2, 3 }, { 1, 3 }, { 1, 2 }, { 0, 3 }, { 0, 2 }, { 0, 1 } };

// UFC: Edge - Non-Incident Vertices
uint const TetrahedronCell::ENV[6][2] =
{ { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 }, { 2, 3 } };

// UFC: Face - Incident Vertices
uint const TetrahedronCell::FIV[4][3] =
{ { 1, 2, 3 }, { 0, 2, 3 }, { 0, 1, 3 }, { 0, 1, 2 } };

// UFC: Face - Non-Incident Vertices
uint const TetrahedronCell::FNV[4][1] =
{ { 0 }, { 1 }, { 2 }, { 3 } };

//-----------------------------------------------------------------------------
TetrahedronCell::TetrahedronCell() :
    CellType("tetrahedron", CellType::tetrahedron, CellType::triangle)
{
}
//-----------------------------------------------------------------------------
void TetrahedronCell::create_entities(uint** e, uint dim, uint const* v) const
{
  // We only need to know how to create edges and faces
  switch (dim)
    {
    case 1:
      // Create the six edges
      e[0][0] = v[2];
      e[0][1] = v[3];
      e[1][0] = v[1];
      e[1][1] = v[3];
      e[2][0] = v[1];
      e[2][1] = v[2];
      e[3][0] = v[0];
      e[3][1] = v[3];
      e[4][0] = v[0];
      e[4][1] = v[2];
      e[5][0] = v[0];
      e[5][1] = v[1];
      break;
    case 2:
      // Create the four faces
      e[0][0] = v[1];
      e[0][1] = v[2];
      e[0][2] = v[3];
      e[1][0] = v[0];
      e[1][1] = v[2];
      e[1][2] = v[3];
      e[2][0] = v[0];
      e[2][1] = v[1];
      e[2][2] = v[3];
      e[3][0] = v[0];
      e[3][1] = v[1];
      e[3][2] = v[2];
      break;
    default:
      error("Invalid topological dimension for creation of entities: %d.", dim);
      break;
    }
}
//-----------------------------------------------------------------------------
void TetrahedronCell::order_entities(MeshTopology& topology, uint i) const
{
  // Sort i - j for i > j: 1 - 0, 2 - 0, 2 - 1, 3 - 0, 3 - 1, 3 - 2
  dolfin_assert(topology.type().cellType() == this->cell_type);

  // Sort local vertices on edges in ascending order, connectivity 1 - 0
  if (topology.connectivity(1, 0))
  {
    dolfin_assert(topology.connectivity(3, 1));

    // Get edges
    Array<uint> const & cell_edges = topology(3, 1)[i];

    // Sort vertices on each edge
    for (uint i = 0; i < 6; ++i)
    {
      Array<uint> & edge_vertices = topology(1, 0)[cell_edges[i]];
      std::sort(edge_vertices.data(), edge_vertices.data() + 2);
    }
  }

  // Sort local vertices on facets in ascending order, connectivity 2 - 0
  if (topology.connectivity(2, 0))
  {
    dolfin_assert(topology.connectivity(3, 2));

    // Get facets
    Array<uint> const & cell_facets = topology(3, 2)[i];

    // Sort vertices on each facet
    for (uint i = 0; i < 4; ++i)
    {
      Array<uint> & facet_vertices = topology(2, 0)[cell_facets[i]];
      std::sort(facet_vertices.data(), facet_vertices.data() + 3);
    }
  }

  // Sort local edges on local facets after non-incident vertex, connectivity 2 - 1
  if (topology.connectivity(2, 1))
  {
    dolfin_assert(topology.connectivity(3, 2));
    dolfin_assert(topology.connectivity(2, 0));
    dolfin_assert(topology.connectivity(1, 0));

    // Get facet numbers
    Array<uint> const & cell_facets = topology(3, 2)[i];

    // Loop over facets on cell
    for (uint i = 0; i < 4; ++i)
    {
      // For each facet number get the global vertex numbers
      Array<uint> const & facet_vertices = topology(2, 0)[cell_facets[i]];

      // For each facet number get the global edge number
      Array<uint> & cell_edges = topology(2, 1)[cell_facets[i]];

      // Loop over vertices on facet
      uint m = 0;
      for (uint j = 0; j < 3; ++j)
      {
        // Loop edges on facet
        for (uint k(m); k < 3; ++k)
        {
          // For each edge number get the global vertex numbers
          Array<uint> const & edge_vertices = topology(1, 0)[cell_edges[k]];

          // Check if the jth vertex of facet i is non-incident on edge k
          if (!std::count(edge_vertices.data(), edge_vertices.data() + 2,
                          facet_vertices[j]))
          {
            // Swap facet numbers
            uint tmp = cell_edges[m];
            cell_edges[m] = cell_edges[k];
            cell_edges[k] = tmp;
            m++;
            break;
          }
        }
      }
    }
  }

  // Sort local vertices on cell in ascending order, connectivity 3 - 0
  if (topology.connectivity(3, 0))
  {
    Array<uint> & cell_vertices = topology(3, 0)[i];
    std::sort(cell_vertices.data(), cell_vertices.data() + 4);
  }

  // Sort local edges on cell after non-incident vertex tuple, connectivity 3 - 1
  if (topology.connectivity(3, 1))
  {
    dolfin_assert(topology.connectivity(1, 0));

    // Get cell vertices and edge numbers
    Array<uint> const & cell_vertices = topology(3, 0)[i];
    Array<uint> & cell_edges = topology(3, 1)[i];

    // Loop two vertices on cell as a lexicographical tuple
    // (i, j): (0,1) (0,2) (0,3) (1,2) (1,3) (2,3)
    uint m = 0;
    for (uint i = 0; i < 3; ++i)
    {
      for (uint j = i + 1; j < 4; ++j)
      {
        // Loop edge numbers
        for (uint k = m; k < 6; ++k)
        {
          // Get local vertices on edge
          Array<uint> const & edge_vertices = topology(1, 0)[cell_edges[k]];

          // Check if the ith and jth vertex of the cell are non-incident on edge k
          if (!std::count(edge_vertices.data(), edge_vertices.data() + 2,
                          cell_vertices[i])
              && !std::count(edge_vertices.data(), edge_vertices.data() + 2,
                             cell_vertices[j]))
          {
            // Swap edge numbers
            uint tmp = cell_edges[m];
            cell_edges[m] = cell_edges[k];
            cell_edges[k] = tmp;
            m++;
            break;
          }
        }
      }
    }
  }

  // Sort local facets on cell after non-incident vertex, connectivity 3 - 2
  if (topology.connectivity(3, 2))
  {
    dolfin_assert(topology.connectivity(2, 0));

    // Get cell vertices and facet numbers
    Array<uint> const & cell_vertices = topology(3, 0)[i];
    Array<uint> & cell_facets = topology(3, 2)[i];

    // Loop vertices on cell
    for (uint i = 0; i < 4; ++i)
    {
      // Loop facets on cell
      for (uint j = i; j < 4; ++j)
      {
        Array<uint> const & facet_vertices = topology(2, 0)[cell_facets[j]];

        // Check if the ith vertex of the cell is non-incident on facet j
        if (!std::count(facet_vertices.data(), facet_vertices.data() + 3,
                        cell_vertices[i]))
        {
          // Swap facet numbers
          uint tmp = cell_facets[i];
          cell_facets[i] = cell_facets[j];
          cell_facets[j] = tmp;
          break;
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
void TetrahedronCell::order_facet(uint vertices[], Facet& facet) const
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
  Point p2 = mesh.geometry().point(facet.entities(0)[2]);
  Point v1 = p1 - p0;
  Point v2 = p2 - p0;
  Point n = v1.cross(v2);

  if (n.dot(p0 - p) < 0.0)
  {
    std::swap(vertices[0], vertices[1]);
  }
}
//-----------------------------------------------------------------------------
void TetrahedronCell::refine_cell(Cell& cell, MeshEditor& editor,
                                 uint& current_cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get vertices and edges
  Array<uint> const & v = cell.entities(0);
  dolfin_assert(!v.empty());
  Array<uint> const & e = cell.entities(1);
  dolfin_assert(!e.empty());

  // Compute indices for the ten new vertices
  uint const v0 = v[0];
  uint const v1 = v[1];
  uint const v2 = v[2];
  uint const v3 = v[3];
  uint const offset = cell.mesh().size(0);
  uint const e0 = offset + e[findEdge(0, cell)];
  uint const e1 = offset + e[findEdge(1, cell)];
  uint const e2 = offset + e[findEdge(2, cell)];
  uint const e3 = offset + e[findEdge(3, cell)];
  uint const e4 = offset + e[findEdge(4, cell)];
  uint const e5 = offset + e[findEdge(5, cell)];

  // Regular refinement: 8 new cells
  uint const cv0[4] = { v0, e3, e4, e5 };
  editor.add_cell(current_cell++, &cv0[0]);
  uint const cv1[4] = { v1, e1, e2, e5 };
  editor.add_cell(current_cell++, &cv1[0]);
  uint const cv2[4] = { v2, e0, e2, e4 };
  editor.add_cell(current_cell++, &cv2[0]);
  uint const cv3[4] = { v3, e0, e1, e3 };
  editor.add_cell(current_cell++, &cv3[0]);
  uint const cv4[4] = { e0, e1, e2, e5 };
  editor.add_cell(current_cell++, &cv4[0]);
  uint const cv5[4] = { e0, e1, e3, e5 };
  editor.add_cell(current_cell++, &cv5[0]);
  uint const cv6[4] = { e0, e2, e4, e5 };
  editor.add_cell(current_cell++, &cv6[0]);
  uint const cv7[4] = { e0, e3, e4, e5 };
  editor.add_cell(current_cell++, &cv7[0]);
}
//-----------------------------------------------------------------------------
auto TetrahedronCell::intersects(MeshEntity const& e, Point const& p) const -> bool
{
  // Adapted from gts_point_is_in_triangle from GTS
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.num_entities(0) == NE[3][0]);

  // Get global index of vertices of the tetrahedron
  uint const ort = orientation((Cell&) e);
  uint const v0 = e.entities(0)[0];
  uint const v1 = e.entities(0)[1];
  uint const v2 = ( ort == 0 ? e.entities(0)[2] : e.entities(0)[3] );
  uint const v3 = ( ort == 0 ? e.entities(0)[3] : e.entities(0)[2] );

  // Get the coordinates of the four vertices
  MeshGeometry const& geometry = e.mesh().geometry();
  real const* x0 = geometry.x(v0);
  real const* x1 = geometry.x(v1);
  real const* x2 = geometry.x(v2);
  real const* x3 = geometry.x(v3);

  // Test orientation of p w.r.t. each face
  real tol = geometry.abs_tolerance(3);
  real d1 = orient3d(x2, x1, x0, &p[0]);
  if (d1 < (-tol)) return false;
  real d2 = orient3d(x0, x3, x2, &p[0]);
  if (d2 < (-tol)) return false;
  real d3 = orient3d(x0, x1, x3, &p[0]);
  if (d3 < (-tol)) return false;
  real d4 = orient3d(x1, x2, x3, &p[0]);
  if (d4 < (-tol)) return false;

  return true;
}
//-----------------------------------------------------------------------------
auto TetrahedronCell::intersects(MeshEntity const& e, Point const&,
                                 Point const&) const -> bool
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.num_entities(0) == NE[3][0]);
  MAYBE_UNUSED(e);

  error("Collision of tetrahedron with segment not implemented");
  return false;
}
//-----------------------------------------------------------------------------
void TetrahedronCell::create_reference_cell(Mesh& mesh) const
{
  MeshEditor me(mesh, CellType::tetrahedron, 3, DOLFIN_COMM_SELF);
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
auto TetrahedronCell::description() const -> std::string
{
  return std::string("tetrahedron (simplex of topological dimension 3)");
}
//-----------------------------------------------------------------------------
void TetrahedronCell::disp() const
{
  section("TetrahedronCell");
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
auto TetrahedronCell::check(Cell& cell) const -> bool
{
  bool ret = CellType::check(cell);

  // UFC convention: cell -> vertices in ascending order
  // These connectivities should always exist, catching assertion if it is not
  // the case is the right behaviour
  Array<uint> const & cell_verts = cell.entities(0);
  dolfin_assert(!cell_verts.empty());
  uint const num_cell_verts = this->num_vertices(this->dim());
  if(!is_sorted(cell_verts.data(), cell_verts.data() + num_cell_verts))
  {
    ret = false;
    warning("CellType : cell vertices are not in ascending order\n"
            "=> cell index = %d", cell.index());
  }

  //
  MeshTopology const& topology = cell.mesh().topology();
  Array<uint> const & v = cell.entities(0);
  dolfin_assert(!v.empty());

  // Check edge -> incident vertices mapping
  if (topology.connectivity(1, 0))
  {
    Array<uint> const & e = cell.entities(1);
    dolfin_assert(!e.empty());
    for (uint i = 0; i < 6; ++i)
    {
      Array<uint> const & ev = topology(1, 0)[e[i]];
      dolfin_assert(!ev.empty());
      for (uint j = 0; j < 2; ++j)
      {
        if (ev[j] != v[EIV[i][j]])
        {
          ret = false;
          warning("CellType : invalid edge -> incident vertices mapping");
        }
      }
    }
  }

  // Check face -> incident vertices mapping
  if (topology.connectivity(2, 0))
  {
    Array<uint> const & f = cell.entities(2);
    dolfin_assert(!f.empty());
    for (uint i = 0; i < 4; ++i)
    {
      Array<uint> const & fv = topology(2, 0)[f[i]];
      dolfin_assert(!fv.empty());
      for (uint j = 0; j < 3; ++j)
      {
        if (fv[j] != v[FIV[i][j]])
        {
          ret = false;
          warning("CellType : invalid face -> incident vertices mapping");
        }
      }
    }
  }

  return ret;
}
//-----------------------------------------------------------------------------
auto TetrahedronCell::findEdge(uint i, Cell const& cell) const -> uint
{
  // Ordering convention for edges (order of non-incident vertices)

  // Get vertices and edges
  Array<uint> const & v = cell.entities(0);
  dolfin_assert(!v.empty());
  Array<uint> const & e = cell.entities(1);
  dolfin_assert(!e.empty());

  // Look for edge satisfying ordering convention
  MeshTopology const& topology = cell.mesh().topology();
  uint const v0 = v[ENV[i][0]];
  uint const v1 = v[ENV[i][1]];
  for (uint j = 0; j < 6; ++j)
  {
    Array<uint> const & ev = topology(1, 0)[e[j]];
    dolfin_assert(!ev.empty());
    if (ev[0] != v0 && ev[0] != v1 && ev[1] != v0 && ev[1] != v1)
    {
      return j;
    }
  }

  // We should not reach this
  error("Unable to find edge with index %d in tetrahedron.", cell.index());

  return 0;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
