// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman 2006.
// Modified by Garth N. Wells 2006.
// Modified by Kristian Oelgaard 2006.
// Modified by Aurelien Larcher, 2015.
//
// First added:  2006-06-05
// Last changed: 2008-06-20

#include <dolfin/mesh/TetrahedronCell.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/GeometricPredicates.h>

#include <algorithm>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
uint const TetrahedronCell::NE[4] =
{ 4, 6, 4, 1 };

// UFC: Number of Vertices (per entity)
uint const TetrahedronCell::NV[4] =
{ 1, 2, 3, 4 };

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
uint const TetrahedronCell::FNV[6][1] =
{ { 0 }, { 1 }, { 2 }, { 3 } };

//-----------------------------------------------------------------------------
TetrahedronCell::TetrahedronCell() :
    CellType(CellType::tetrahedron, CellType::triangle)
{
}
//-----------------------------------------------------------------------------
TetrahedronCell::~TetrahedronCell()
{
}
//-----------------------------------------------------------------------------
uint TetrahedronCell::dim() const
{
  return 3;
}
//-----------------------------------------------------------------------------
uint TetrahedronCell::numEntities(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NE[dim];
}
//-----------------------------------------------------------------------------
uint TetrahedronCell::numVertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NV[dim];
}
//-----------------------------------------------------------------------------
uint TetrahedronCell::orientation(Cell const& cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = cell.mesh().geometry();
  uint const * vertices = cell.entities(0);
  real const * v0 = geometry.x(vertices[0]);
  real const * v1 = geometry.x(vertices[1]);
  real const * v2 = geometry.x(vertices[2]);
  real const * v3 = geometry.x(vertices[3]);

  real a = +((v1[1] - v0[1]) * (v2[2] - v0[2])
      - (v1[2] - v0[2]) * (v2[1] - v0[1])) * (v3[0] - v0[0])
      + ((v1[2] - v0[2]) * (v2[0] - v0[0]) - (v1[0] - v0[0]) * (v2[2] - v0[2]))
          * (v3[1] - v0[1])
      + ((v1[0] - v0[0]) * (v2[1] - v0[1]) - (v1[1] - v0[1]) * (v2[0] - v0[0]))
          * (v3[2] - v0[2]);

  return (a < 0.0 ? 1 : 0);
}
//-----------------------------------------------------------------------------
void TetrahedronCell::createEntities(uint** e, uint dim, uint const* v) const
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
void TetrahedronCell::orderEntities(Cell& cell) const
{
  // Sort i - j for i > j: 1 - 0, 2 - 0, 2 - 1, 3 - 0, 3 - 1, 3 - 2
  dolfin_assert(cell.type() == this->cell_type);

  // Get mesh topology
  MeshTopology& topology = cell.mesh().topology();

  // Sort local vertices on edges in ascending order, connectivity 1 - 0
  if (topology(1, 0).size() > 0)
  {
    dolfin_assert(topology(3, 1).size() > 0);

    // Get edges
    uint* cell_edges = cell.entities(1);

    // Sort vertices on each edge
    for (uint i = 0; i < 6; ++i)
    {
      uint* edge_vertices = topology(1, 0)(cell_edges[i]);
      std::sort(edge_vertices, edge_vertices + 2);
    }
  }

  // Sort local vertices on facets in ascending order, connectivity 2 - 0
  if (topology(2, 0).size() > 0)
  {
    dolfin_assert(topology(3, 2).size() > 0);

    // Get facets
    uint* cell_facets = cell.entities(2);

    // Sort vertices on each facet
    for (uint i = 0; i < 4; ++i)
    {
      uint* facet_vertices = topology(2, 0)(cell_facets[i]);
      std::sort(facet_vertices, facet_vertices + 3);
    }
  }

  // Sort local edges on local facets after non-incident vertex, connectivity 2 - 1
  if (topology(2, 1).size() > 0)
  {
    dolfin_assert(topology(3, 2).size() > 0);
    dolfin_assert(topology(2, 0).size() > 0);
    dolfin_assert(topology(1, 0).size() > 0);

    // Get facet numbers
    uint* cell_facets = cell.entities(2);

    // Loop over facets on cell
    for (uint i = 0; i < 4; ++i)
    {
      // For each facet number get the global vertex numbers
      uint* facet_vertices = topology(2, 0)(cell_facets[i]);

      // For each facet number get the global edge number
      uint* cell_edges = topology(2, 1)(cell_facets[i]);

      // Loop over vertices on facet
      uint m = 0;
      for (uint j = 0; j < 3; ++j)
      {
        // Loop edges on facet
        for (uint k(m); k < 3; ++k)
        {
          // For each edge number get the global vertex numbers
          uint* edge_vertices = topology(1, 0)(cell_edges[k]);

          // Check if the jth vertex of facet i is non-incident on edge k
#if __SUNPRO_CC
          int n1 = 0;
          std::count(edge_vertices, edge_vertices+2, facet_vertices[j], n1);
          if (!n1)
#else
          if (!std::count(edge_vertices, edge_vertices + 2, facet_vertices[j]))
#endif
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
  if (topology(3, 0).size() > 0)
  {
    uint* cell_vertices = cell.entities(0);
    std::sort(cell_vertices, cell_vertices + 4);
  }

  // Sort local edges on cell after non-incident vertex tuple, connectivity 3 - 1
  if (topology(3, 1).size() > 0)
  {
    dolfin_assert(topology(1, 0).size() > 0);

    // Get cell vertices and edge numbers
    uint* cell_vertices = cell.entities(0);
    uint* cell_edges = cell.entities(1);

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
          uint* edge_vertices = topology(1, 0)(cell_edges[k]);

          // Check if the ith and jth vertex of the cell are non-incident on edge k
#if __SUNPRO_CC
          int n1 = 0;
          int n2 = 0;
          std::count(edge_vertices, edge_vertices+2, cell_vertices[i], n1);
          std::count(edge_vertices, edge_vertices+2, cell_vertices[j], n2);
          if (!n1 && !n2 )
#else
          if (!std::count(edge_vertices, edge_vertices + 2, cell_vertices[i])
              && !std::count(edge_vertices, edge_vertices + 2,
                             cell_vertices[j]))
#endif
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
  if (topology(3, 2).size() > 0)
  {
    dolfin_assert(topology(2, 0).size() > 0);

    // Get cell vertices and facet numbers
    uint* cell_vertices = cell.entities(0);
    uint* cell_facets = cell.entities(2);

    // Loop vertices on cell
    for (uint i = 0; i < 4; ++i)
    {
      // Loop facets on cell
      for (uint j = i; j < 4; ++j)
      {
        uint* facet_vertices = topology(2, 0)(cell_facets[j]);

        // Check if the ith vertex of the cell is non-incident on facet j
#if __SUNPRO_CC
        int n1 = 0;
        std::count(facet_vertices, facet_vertices+3, cell_vertices[i], n1);
        if (!n1)
#else
        if (!std::count(facet_vertices, facet_vertices + 3, cell_vertices[i]))
#endif
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
void TetrahedronCell::refine_cell(Cell& cell, MeshEditor& editor,
                                 uint& current_cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get vertices and edges
  uint const* v = cell.entities(0);
  dolfin_assert(v);
  uint const* e = cell.entities(1);
  dolfin_assert(e);

  // Compute indices for the ten new vertices
  uint const v0 = v[0];
  uint const v1 = v[1];
  uint const v2 = v[2];
  uint const v3 = v[3];
  uint const offset = cell.mesh().numVertices();
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
uint TetrahedronCell::num_refined_cells() const
{
  return 8;
}
//-----------------------------------------------------------------------------
uint TetrahedronCell::num_refined_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return (dim > 1 ? 0 : 1);
}
//-----------------------------------------------------------------------------
bool TetrahedronCell::refinement_needs_entities(uint dim) const
{
  dolfin_assert(dim <= TD);
  if (dim > 1)
  {
    return false;
  }
  return true;
}
//-----------------------------------------------------------------------------
real TetrahedronCell::volume(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.numEntities(0) == NE[0]);

  // Get the coordinates of the four vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const* vertices = entity.entities(0);
  real const* x0 = geometry.x(vertices[0]);
  real const* x1 = geometry.x(vertices[1]);
  real const* x2 = geometry.x(vertices[2]);
  real const* x3 = geometry.x(vertices[3]);

  // Formula for volume from http://mathworld.wolfram.com
  real V = (+x0[0]
      * (x1[1] * x2[2] + x3[1] * x1[2] + x2[1] * x3[2] - x2[1] * x1[2]
          - x1[1] * x3[2] - x3[1] * x2[2])
      - x1[0]
          * (x0[1] * x2[2] + x3[1] * x0[2] + x2[1] * x3[2] - x2[1] * x0[2]
              - x0[1] * x3[2] - x3[1] * x2[2])
      + x2[0]
          * (x0[1] * x1[2] + x3[1] * x0[2] + x1[1] * x3[2] - x1[1] * x0[2]
              - x0[1] * x3[2] - x3[1] * x1[2])
      - x3[0]
          * (x0[1] * x1[2] + x1[1] * x2[2] + x2[1] * x0[2] - x1[1] * x0[2]
              - x2[1] * x1[2] - x0[1] * x2[2]));

  return std::abs(V) / 6.0;
}
//-----------------------------------------------------------------------------
real TetrahedronCell::diameter(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.numEntities(0) == NE[0]);

  // Get the coordinates of the four vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const* vertices = entity.entities(0);
  real const* x0 = geometry.x(vertices[0]);
  real const* x1 = geometry.x(vertices[1]);
  real const* x2 = geometry.x(vertices[2]);
  real const* x3 = geometry.x(vertices[3]);

  // Compute edge lengths
  real a = 0.0;
  real b = 0.0;
  real c = 0.0;
  real aa = 0.0;
  real bb = 0.0;
  real cc = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    a += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    b += (x0[i] - x2[i]) * (x0[i] - x2[i]);
    c += (x0[i] - x1[i]) * (x0[i] - x1[i]);
    aa += (x0[i] - x3[i]) * (x0[i] - x3[i]);
    bb += (x1[i] - x3[i]) * (x1[i] - x3[i]);
    cc += (x2[i] - x3[i]) * (x2[i] - x3[i]);
  }

  real hmax = a;
  hmax = std::max(b, hmax);
  hmax = std::max(c, hmax);
  hmax = std::max(aa, hmax);
  hmax = std::max(bb, hmax);
  hmax = std::max(cc, hmax);
  return std::sqrt(hmax);
}
//-----------------------------------------------------------------------------
real TetrahedronCell::circumradius(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.numEntities(0) == NE[0]);

  // Get the coordinates of the four vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const* vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);

  // Compute edge lengths
  real a = 0.0;
  real b = 0.0;
  real c = 0.0;
  real aa = 0.0;
  real bb = 0.0;
  real cc = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    a += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    b += (x0[i] - x2[i]) * (x0[i] - x2[i]);
    c += (x0[i] - x1[i]) * (x0[i] - x1[i]);
    aa += (x0[i] - x3[i]) * (x0[i] - x3[i]);
    bb += (x1[i] - x3[i]) * (x1[i] - x3[i]);
    cc += (x2[i] - x3[i]) * (x2[i] - x3[i]);
  }

  // Compute "area" of triangle with strange side lengths
  real la = a * aa;
  real lb = b * bb;
  real lc = c * cc;
  real s = 0.5 * (la + lb + lc);
  real area = std::sqrt(s * (s - la) * (s - lb) * (s - lc));

  // Formula for volume from http://mathworld.wolfram.com
  real V = (+x0[0]
      * (x1[1] * x2[2] + x3[1] * x1[2] + x2[1] * x3[2] - x2[1] * x1[2]
          - x1[1] * x3[2] - x3[1] * x2[2])
      - x1[0]
          * (x0[1] * x2[2] + x3[1] * x0[2] + x2[1] * x3[2] - x2[1] * x0[2]
              - x0[1] * x3[2] - x3[1] * x2[2])
      + x2[0]
          * (x0[1] * x1[2] + x3[1] * x0[2] + x1[1] * x3[2] - x1[1] * x0[2]
              - x0[1] * x3[2] - x3[1] * x1[2])
      - x3[0]
          * (x0[1] * x1[2] + x1[1] * x2[2] + x2[1] * x0[2] - x1[1] * x0[2]
              - x2[1] * x1[2] - x0[1] * x2[2]));

  // Formula for circumradius from http://mathworld.wolfram.com
  return area / (6.0 * V);
}
//-----------------------------------------------------------------------------
Point TetrahedronCell::midpoint(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.numEntities(0) == NE[0]);

  // Get the coordinates of the vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const* vertices = entity.entities(0);
  real const* x0 = geometry.x(vertices[0]);
  real const* x1 = geometry.x(vertices[1]);
  real const* x2 = geometry.x(vertices[2]);
  real const* x3 = geometry.x(vertices[3]);
  Point p;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    p[i] = 0.25 * ( x0[i] + x1[i] + x2[i] + x3[i] );
  }
  return p;
}
//-----------------------------------------------------------------------------
Point TetrahedronCell::normal(Cell const& cell, uint facet) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Create facet from the mesh and local facet number
  Cell& c = const_cast<Cell&>(cell);
  Facet f(c.mesh(), c.entities(2)[facet]);

  // Get global index of opposite vertex
  uint const v0 = cell.entities(0)[facet];

  // Get global index of vertices on the facet
  uint v1 = f.entities(0)[0];
  uint v2 = f.entities(0)[1];
  uint v3 = f.entities(0)[2];

  // Get the coordinates of the four vertices
  MeshGeometry const& geometry = cell.mesh().geometry();
  real const* x0 = geometry.x(v0);
  real const* x1 = geometry.x(v1);
  real const* x2 = geometry.x(v2);
  real const* x3 = geometry.x(v3);

  // Create vectors
  Point e0;
  Point e1;
  Point e2;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    e0[i] = (x0[i] - x1[i]);
    e1[i] = (x2[i] - x1[i]);
    e2[i] = (x3[i] - x1[i]);
  }

  // Compute normal vector
  Point n = e1.cross(e2);

  // Normalize
  n /= n.norm();

  // Flip direction of normal so it points outward
  if (n.dot(e0) > 0) n *= -1.0;

  return n;
}
//-----------------------------------------------------------------------------
dolfin::real TetrahedronCell::facetArea(Cell const& cell, uint facet) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Create facet from the mesh and local facet number
  Cell& c = const_cast<Cell&>(cell);
  Facet f(c.mesh(), c.entities(2)[facet]);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = cell.mesh().geometry();
  uint const* vertices = cell.entities(0);
  real const* x0 = geometry.x(vertices[0]);
  real const* x1 = geometry.x(vertices[1]);
  real const* x2 = geometry.x(vertices[2]);

  // Compute area of triangle embedded in R^3
  real v0 = (x0[1] * x1[2] + x0[2] * x2[1] + x1[1] * x2[2])
      - (x2[1] * x1[2] + x2[2] * x0[1] + x1[1] * x0[2]);
  real v1 = (x0[2] * x1[0] + x0[0] * x2[2] + x1[2] * x2[0])
      - (x2[2] * x1[0] + x2[0] * x0[2] + x1[2] * x0[0]);
  real v2 = (x0[0] * x1[1] + x0[1] * x2[0] + x1[0] * x2[1])
      - (x2[0] * x1[1] + x2[1] * x0[0] + x1[0] * x0[1]);

  // Formula for area from http://mathworld.wolfram.com
  return 0.5 * std::sqrt(v0 * v0 + v1 * v1 + v2 * v2);
}
//-----------------------------------------------------------------------------
bool TetrahedronCell::intersects(MeshEntity const& e, Point const& p) const
{
  // Adapted from gts_point_is_in_triangle from GTS
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.numEntities(0) == NE[0]);

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
bool TetrahedronCell::intersects(MeshEntity const& e, Point const& p1,
                                 Point const& p2) const
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.numEntities(0) == NE[0]);

  error("Collision of tetrahedron with segment not implemented");
  return false;
}
//-----------------------------------------------------------------------------
Mesh TetrahedronCell::create_reference_cell() const
{
  Mesh refcell;
  MeshEditor me(refcell, CellType::tetrahedron, 3);
  me.init_vertices(4);
  me.add_vertex(0, VC[0]);
  me.add_vertex(1, VC[1]);
  me.add_vertex(2, VC[2]);
  me.add_vertex(3, VC[3]);
  me.init_cells(1);
  uint const cv0[4] = { 0, 1, 2, 3 };
  me.add_cell(0, cv0);
  me.close();
  return refcell;
}
//-----------------------------------------------------------------------------
std::string TetrahedronCell::description() const
{
  return std::string("tetrahedron (simplex of topological dimension 3)");
}
//-----------------------------------------------------------------------------
void TetrahedronCell::disp() const
{
  message("TetrahedronCell");
  begin(  "---------------");
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
void TetrahedronCell::check(Cell& cell) const
{
  CellType::check(cell);

  //
  MeshTopology const& topology = cell.mesh().topology();
  uint const* v = cell.entities(0);
  dolfin_assert(v);

  // Check edge -> incident vertices mapping
  uint const* e = cell.entities(1);
  dolfin_assert(e);
  for (uint i = 0; i < 6; ++i)
  {
    uint const * ev = topology(1, 0)(e[i]);
    dolfin_assert(ev);
    for (uint j = 0; j < 2; ++j)
    {
      if (ev[j] != v[EIV[i][j]])
      {
        error("CellType::check : invalid edge -> incident vertices mapping");
      }
    }
  }

  // Check face -> incident vertices mapping
  uint const* f = cell.entities(2);
  dolfin_assert(f);
  for (uint i = 0; i < 4; ++i)
  {
    uint const * fv = topology(2, 0)(f[i]);
    dolfin_assert(fv);
    for (uint j = 0; j < 3; ++j)
    {
      if (fv[j] != v[FIV[i][j]])
      {
        error("CellType::check : invalid face -> incident vertices mapping");
      }
    }
  }
}
//-----------------------------------------------------------------------------
uint TetrahedronCell::findEdge(uint i, Cell const& cell) const
{
  // Ordering convention for edges (order of non-incident vertices)

  // Get vertices and edges
  uint const* v = cell.entities(0);
  dolfin_assert(v);
  uint const* e = cell.entities(1);
  dolfin_assert(e);

  // Look for edge satisfying ordering convention
  MeshTopology const& topology = cell.mesh().topology();
  uint const v0 = v[ENV[i][0]];
  uint const v1 = v[ENV[i][1]];
  for (uint j = 0; j < 6; ++j)
  {
    uint const* ev = topology(1, 0)(e[j]);
    dolfin_assert(ev);
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
