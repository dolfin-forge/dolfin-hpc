// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
//

#include <dolfin/mesh/HexahedronCell.h>

#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Cell.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------

// UFC: Number of Entities
uint const HexahedronCell::NE[4] =
{ 8, 12, 6, 1 };

// UFC: Number of Vertices (per entity)
uint const HexahedronCell::NV[4] =
{ 1, 2, 4, 8 };

// UFC: Vertex Coordinates
real const HexahedronCell::VC[8][3] =
{ { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 },
  { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 } };

// UFC: Edge - Incident Vertices
uint const HexahedronCell::EIV[12][2] =
{ { 6, 7 }, { 5, 6 }, { 4, 7 }, { 4, 5 }, { 3, 7 }, { 2, 6 },
  { 2, 3 }, { 1, 5 }, { 1, 2 }, { 0, 4 }, { 0, 3 }, { 0, 1 } };

// UFC: Edge - Non-Incident Vertices
uint const HexahedronCell::ENV[12][6] =
{ { 0, 1, 2, 3, 4, 5 }, { 0, 1, 2, 3, 4, 7 }, { 0, 1, 2, 3, 5, 6 },
  { 0, 1, 2, 3, 6, 7 }, { 0, 1, 2, 4, 5, 6 }, { 0, 1, 3, 4, 5, 7 },
  { 0, 1, 4, 5, 6, 7 }, { 0, 2, 3, 4, 6, 7 }, { 0, 3, 4, 5, 6, 7 },
  { 1, 2, 3, 5, 6, 7 }, { 1, 2, 4, 5, 6, 7 }, { 2, 3, 4, 5, 6, 7 } };

// UFC: Face - Incident Vertices
uint const HexahedronCell::FIV[6][4] =
{ { 4, 5, 6, 7 }, { 2, 3, 6, 7 }, { 1, 2, 5, 6 },
  { 0, 3, 4, 7 }, { 0, 1, 4, 5 }, { 0, 1, 2, 3 } };

// UFC: Face - Non-Incident Vertices
uint const HexahedronCell::FNV[6][4] =
{ { 0, 1, 2, 3 }, { 0, 1, 4, 5 }, { 0, 3, 4, 7 },
  { 1, 2, 5, 6 }, { 2, 3, 6, 7 }, { 4, 5, 6, 7 } };

//-----------------------------------------------------------------------------
HexahedronCell::HexahedronCell() :
    CellType(CellType::hexahedron, CellType::quadrilateral)
{
}
//-----------------------------------------------------------------------------
HexahedronCell::~HexahedronCell()
{
}
//-----------------------------------------------------------------------------
uint HexahedronCell::dim() const
{
  return 3;
}
//-----------------------------------------------------------------------------
uint HexahedronCell::num_entities(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NE[dim];
}
//-----------------------------------------------------------------------------
uint HexahedronCell::num_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NV[dim];
}
//-----------------------------------------------------------------------------
uint HexahedronCell::orientation(Cell const& cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get the coordinates of vertices v0, v1, v3 and v4
  MeshGeometry const& geometry = cell.mesh().geometry();
  uint const * vertices = cell.entities(0);
  real const * v0 = geometry.x(vertices[0]);
  real const * v1 = geometry.x(vertices[1]);
  real const * v3 = geometry.x(vertices[3]);
  real const * v4 = geometry.x(vertices[4]);

  // Check whether (v0v1, v0v3, v0v4) is counter-clockwise
  real const S =
    + ((v1[1] - v0[1]) * (v3[2] - v0[2]) - (v1[2] - v0[2]) * (v3[1] - v0[1]))
          * (v4[0] - v0[0])
    + ((v1[2] - v0[2]) * (v3[0] - v0[0]) - (v1[0] - v0[0]) * (v3[2] - v0[2]))
          * (v4[1] - v0[1])
    + ((v1[0] - v0[0]) * (v3[1] - v0[1]) - (v1[1] - v0[1]) * (v3[0] - v0[0]))
          * (v4[2] - v0[2]);
  return (S < 0.0 ? 1 : 0);
}
//-----------------------------------------------------------------------------
void HexahedronCell::create_entities(uint** e, uint dim, uint const* v) const
{
  // We only need to know how to create edges and faces
  switch (dim)
    {
    case 1:
      // Create the twelve edges
      e[0][0] = v[6];
      e[0][1] = v[7];
      e[1][0] = v[5];
      e[1][1] = v[6];
      e[2][0] = v[4];
      e[2][1] = v[7];
      e[3][0] = v[4];
      e[3][1] = v[5];
      e[4][0] = v[3];
      e[4][1] = v[7];
      e[5][0] = v[2];
      e[5][1] = v[6];
      e[6][0] = v[2];
      e[6][1] = v[3];
      e[7][0] = v[1];
      e[7][1] = v[5];
      e[8][0] = v[1];
      e[8][1] = v[2];
      e[9][0] = v[0];
      e[9][1] = v[4];
      e[10][0] = v[0];
      e[10][1] = v[3];
      e[11][0] = v[0];
      e[11][1] = v[1];
      break;
    case 2:
      // Create the six faces
      e[0][0] = v[4];
      e[0][1] = v[5];
      e[0][2] = v[6];
      e[0][3] = v[7];
      e[1][0] = v[2];
      e[1][1] = v[3];
      e[1][2] = v[6];
      e[1][3] = v[7];
      e[2][0] = v[1];
      e[2][1] = v[2];
      e[2][2] = v[5];
      e[2][3] = v[6];
      e[3][0] = v[0];
      e[3][1] = v[3];
      e[3][2] = v[4];
      e[3][3] = v[7];
      e[4][0] = v[0];
      e[4][1] = v[1];
      e[4][2] = v[4];
      e[4][3] = v[5];
      e[5][0] = v[0];
      e[5][1] = v[1];
      e[5][2] = v[2];
      e[5][3] = v[3];
      break;
    default:
      error("Invalid topological dimension for creation of entities: %d.", dim);
      break;
    }
}
//-----------------------------------------------------------------------------
void HexahedronCell::order_entities(Cell& cell) const
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
    for (uint i = 0; i < 12; ++i)
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
    for (uint i = 0; i < 6; ++i)
    {
      uint* facet_vertices = topology(2, 0)(cell_facets[i]);
      std::sort(facet_vertices, facet_vertices + 4);
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
    for (uint i = 0; i < 6; ++i)
    {
      // For each facet number get the global vertex numbers
      uint* facet_vertices = topology(2, 0)(cell_facets[i]);

      // For each facet number get the global edge number
      uint* cell_edges = topology(2, 1)(cell_facets[i]);

      // Loop over vertices on facet
      uint m = 0;
      for (uint j = 0; j < 4; ++j)
      {
        // Loop edges on facet
        for (uint k(m); k < 4; ++k)
        {
          // For each edge number get the global vertex numbers
          uint* edge_vertices = topology(1, 0)(cell_edges[k]);

          // Check if the jth vertex of facet i is non-incident on edge k
#if __SUNPRO_CC
          int n1 = 0;
          std::count(edge_vertices, edge_vertices + 2, facet_vertices[j], n1);
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
    std::sort(cell_vertices, cell_vertices + 8);
  }

  // Sort local edges on cell after non-incident vertex tuple, connectivity 3-1
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
          std::count(edge_vertices, edge_vertices + 2, cell_vertices[i], n1);
          std::count(edge_vertices, edge_vertices + 2, cell_vertices[j], n2);
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
        std::count(facet_vertices, facet_vertices + 4, cell_vertices[i], n1);
        if (!n1)
#else
        if (!std::count(facet_vertices, facet_vertices + 4, cell_vertices[i]))
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
void HexahedronCell::refine_cell(Cell& cell, MeshEditor& editor,
                                uint& current_cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get vertices and edges
  uint const * v = cell.entities(0);
  dolfin_assert(v);
  uint const * e = cell.entities(1);
  dolfin_assert(e);
  uint const * f = cell.entities(2);
  dolfin_assert(f);

  // Compute indices for the twenty-seven new vertices
  uint const v00 = v[0];
  uint const v01 = v[1];
  uint const v02 = v[2];
  uint const v03 = v[3];
  uint const v04 = v[4];
  uint const v05 = v[5];
  uint const v06 = v[6];
  uint const v07 = v[7];
  uint const eoffset = cell.mesh().numVertices();
  uint const e00 = eoffset + e[findEdge(0, cell)];
  uint const e01 = eoffset + e[findEdge(1, cell)];
  uint const e02 = eoffset + e[findEdge(2, cell)];
  uint const e03 = eoffset + e[findEdge(3, cell)];
  uint const e04 = eoffset + e[findEdge(4, cell)];
  uint const e05 = eoffset + e[findEdge(5, cell)];
  uint const e06 = eoffset + e[findEdge(6, cell)];
  uint const e07 = eoffset + e[findEdge(7, cell)];
  uint const e08 = eoffset + e[findEdge(8, cell)];
  uint const e09 = eoffset + e[findEdge(9, cell)];
  uint const e10 = eoffset + e[findEdge(10, cell)];
  uint const e11 = eoffset + e[findEdge(11, cell)];
  uint const foffset = eoffset + cell.mesh().numEdges();
  uint const f00 = foffset + f[findFace(0, cell)];
  uint const f01 = foffset + f[findFace(1, cell)];
  uint const f02 = foffset + f[findFace(2, cell)];
  uint const f03 = foffset + f[findFace(3, cell)];
  uint const f04 = foffset + f[findFace(4, cell)];
  uint const f05 = foffset + f[findFace(5, cell)];
  uint const coffset = foffset + cell.mesh().numFaces();
  uint const c00 = coffset + cell.index();

  // Add the eight new cells
  uint const cv0[8] = { v00, e11, f05, e10, e09, f04, c00, f03 };
  editor.add_cell(current_cell++, &cv0[0]);
  uint const cv1[8] = { v01, e08, f05, e11, e07, f02, c00, f04 };
  editor.add_cell(current_cell++, &cv1[0]);
  uint const cv2[8] = { v02, e06, f05, e08, e05, f01, c00, f02 };
  editor.add_cell(current_cell++, &cv2[0]);
  uint const cv3[8] = { v03, e10, f05, e06, e04, f03, c00, f01 };
  editor.add_cell(current_cell++, &cv3[0]);
  uint const cv4[8] = { v04, e03, f00, e02, e09, f04, c00, f03 };
  editor.add_cell(current_cell++, &cv4[0]);
  uint const cv5[8] = { v05, e01, f00, e03, e07, f02, c00, f04 };
  editor.add_cell(current_cell++, &cv5[0]);
  uint const cv6[8] = { v06, e00, f00, e01, e05, f01, c00, f02 };
  editor.add_cell(current_cell++, &cv6[0]);
  uint const cv7[8] = { v07, e02, f00, e00, e04, f03, c00, f01 };
  editor.add_cell(current_cell++, &cv7[0]);
}
//-----------------------------------------------------------------------------
uint HexahedronCell::num_refined_cells() const
{
  return 8;
}
//-----------------------------------------------------------------------------
uint HexahedronCell::num_refined_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return 1;
}
//-----------------------------------------------------------------------------
bool HexahedronCell::refinement_needs_entities(uint dim) const
{
  dolfin_assert(dim <= TD);
  return true;
}
//-----------------------------------------------------------------------------
real HexahedronCell::volume(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.numEntities(0) == NE[0]);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * v0 = geometry.x(vertices[0]);
  real const * v1 = geometry.x(vertices[1]);
  real const * v3 = geometry.x(vertices[3]);
  real const * v4 = geometry.x(vertices[4]);

  // Compute mixed product of (v0v1, v0v3, v0v4)
  real const S =
      + ((v1[1] - v0[1]) * (v3[2] - v0[2]) - (v1[2] - v0[2]) * (v3[1] - v0[1]))
            * (v4[0] - v0[0])
      + ((v1[2] - v0[2]) * (v3[0] - v0[0]) - (v1[0] - v0[0]) * (v3[2] - v0[2]))
            * (v4[1] - v0[1])
      + ((v1[0] - v0[0]) * (v3[1] - v0[1]) - (v1[1] - v0[1]) * (v3[0] - v0[0]))
            * (v4[2] - v0[2]);
  return std::fabs(S);
}
//-----------------------------------------------------------------------------
real HexahedronCell::diameter(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.numEntities(0) == NE[0]);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
  real const * x4 = geometry.x(vertices[4]);
  real const * x5 = geometry.x(vertices[5]);
  real const * x6 = geometry.x(vertices[6]);
  real const * x7 = geometry.x(vertices[7]);

  // Compute maximum diagonal
  real d0 = 0.0;
  real d1 = 0.0;
  real d2 = 0.0;
  real d3 = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    d0 += (x6[i] - x0[i]) * (x6[i] - x0[i]);
    d1 += (x7[i] - x1[i]) * (x7[i] - x1[i]);
    d2 += (x4[i] - x2[i]) * (x4[i] - x2[i]);
    d3 += (x5[i] - x3[i]) * (x5[i] - x3[i]);
  }
  return std::sqrt(std::max(std::max(d0, d1),std::max(d2, d3)));
}
//-----------------------------------------------------------------------------
real HexahedronCell::circumradius(MeshEntity const& entity) const
{
  return this->diameter(entity);
}
//-----------------------------------------------------------------------------
Point HexahedronCell::midpoint(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.numEntities(0) == NE[0]);

  // Get the coordinates of the vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
  real const * x4 = geometry.x(vertices[4]);
  real const * x5 = geometry.x(vertices[5]);
  real const * x6 = geometry.x(vertices[6]);
  real const * x7 = geometry.x(vertices[7]);
  Point p;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    p[i] = 0.125
        * (x0[i] + x1[i] + x2[i] + x3[i] + x4[i] + x5[i] + x6[i] + x7[i]);
  }
  return p;
}
//-----------------------------------------------------------------------------
Point HexahedronCell::normal(Cell const& cell, uint facet) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Create facet from the mesh and local facet number
  Cell& c = const_cast<Cell&>(cell);
  Facet f(c.mesh(), c.entities(1)[facet]);
  MeshGeometry const& geometry = cell.mesh().geometry();

  // first non-incident vertex
  real const * vf = geometry.x(c.entities(0)[FNV[facet][0]]);

  // vertices on the facet
  real const * v0 = geometry.x(f.entities(0)[FIV[facet][0]]);
  real const * v1 = geometry.x(f.entities(0)[FIV[facet][1]]);
  real const * v3 = geometry.x(f.entities(0)[FIV[facet][3]]);

  // Vector normal to facet
  real nx = + (v1[1] - v0[1]) * (v3[2] - v0[2])
            - (v1[2] - v0[2]) * (v3[1] - v0[1]);
  real ny = + (v1[2] - v0[2]) * (v3[0] - v0[0])
            - (v1[0] - v0[0]) * (v3[2] - v0[2]);
  real nz = + (v1[0] - v0[0]) * (v3[1] - v0[1])
            - (v1[1] - v0[1]) * (v3[0] - v0[0]);
  real const nn = std::sqrt(nx * nx + ny * ny + nz * nz);
  nx /= nn;
  ny /= nn;
  nz /= nn;
  Point n(nx, ny, nz);

  // Flip direction of normal so it points outward
  real ps = + n[0] * (vf[0] - v0[0])
            + n[1] * (vf[1] - v0[1])
            + n[2] * (vf[2] - v0[2]);
  if (ps > 0)
  {
    n *= -1.0;
  }
  return n;
}
//-----------------------------------------------------------------------------
real HexahedronCell::facet_area(Cell const& cell, uint facet) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Create facet from the mesh and local facet number
  Cell& c = const_cast<Cell&>(cell);
  Facet f(c.mesh(), c.entities(1)[facet]);

  // Get global index of vertices on the facet
  uint const v0 = f.entities(0)[FIV[facet][0]];
  uint const v1 = f.entities(0)[FIV[facet][1]];
  uint const v2 = f.entities(0)[FIV[facet][2]];
  uint const v3 = f.entities(0)[FIV[facet][3]];

  // Get the coordinates of the two vertices
  MeshGeometry const& geometry = cell.mesh().geometry();
  real const * p = geometry.x(v0);
  real const * q = geometry.x(v1);
  real const * r = geometry.x(v0);
  real const * s = geometry.x(v1);

  // Compute the area with diagonal formula
  real c0 = (r[1] - p[1]) * (s[2] - q[2]) - (r[2] - p[2]) * (s[1] - q[1]);
  real c1 = (r[2] - p[2]) * (s[0] - q[0]) - (r[0] - p[0]) * (s[2] - q[2]);
  real c2 = (r[0] - p[0]) * (s[1] - q[1]) - (r[1] - p[1]) * (s[0] - q[0]);
  return 0.5 * std::sqrt(c0 * c0 + c1 * c1 + c2 * c2);
}
//-----------------------------------------------------------------------------
bool HexahedronCell::intersects(MeshEntity const& e, Point const& p) const
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.numEntities(0) == NE[0]);

  // Get the coordinates of the vertices
  MeshGeometry const& geometry = e.mesh().geometry();
  uint const* vertices = e.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
  real const * x4 = geometry.x(vertices[4]);
  real const * x5 = geometry.x(vertices[5]);
  real const * x6 = geometry.x(vertices[6]);
  real const * x7 = geometry.x(vertices[7]);

  error("Collision of hexahedron with point not implemented.");

  return true;
}
//-----------------------------------------------------------------------------
bool HexahedronCell::intersects(MeshEntity const& e, Point const& p1,
                                Point const& p2) const
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.numEntities(0) == NE[0]);

  // Get the coordinates of the vertices
  MeshGeometry const& geometry = e.mesh().geometry();
  uint const* vertices = e.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real const * x3 = geometry.x(vertices[3]);
  real const * x4 = geometry.x(vertices[4]);
  real const * x5 = geometry.x(vertices[5]);
  real const * x6 = geometry.x(vertices[6]);
  real const * x7 = geometry.x(vertices[7]);

  error("Collision of hexahedron with segment not implemented.");

  return true;
}
//-----------------------------------------------------------------------------
Mesh HexahedronCell::create_reference_cell() const
{
  Mesh refcell;
  MeshEditor me(refcell, CellType::hexahedron, 3);
  me.init_vertices(8);
  me.add_vertex(0, VC[0]);
  me.add_vertex(1, VC[1]);
  me.add_vertex(2, VC[2]);
  me.add_vertex(3, VC[3]);
  me.add_vertex(4, VC[4]);
  me.add_vertex(5, VC[5]);
  me.add_vertex(6, VC[6]);
  me.add_vertex(7, VC[7]);
  me.init_cells(1);
  uint const cv0[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  me.add_cell(0, cv0);
  me.close();
  return refcell;
}
//-----------------------------------------------------------------------------
std::string HexahedronCell::description() const
{
  return std::string("hexahedron (hypercube of topological dimension 3)");
}
//-----------------------------------------------------------------------------
void HexahedronCell::disp() const
{
  message("HexahedronCell");
  begin(  "--------------");
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
void HexahedronCell::check(Cell& cell) const
{
  CellType::check(cell);

  //
  MeshTopology const& topology = cell.mesh().topology();
  uint const* v = cell.entities(0);
  dolfin_assert(v);

  // Check edge -> incident vertices mapping
  uint const* e = cell.entities(1);
  dolfin_assert(e);
  for (uint i = 0; i < 12; ++i)
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
  for (uint i = 0; i < 6; ++i)
  {
    uint const * fv = topology(2, 0)(f[i]);
    dolfin_assert(fv);
    for (uint j = 0; j < 4; ++j)
    {
      if (fv[j] != v[FIV[i][j]])
      {
        error("CellType::check : invalid face -> incident vertices mapping");
      }
    }
  }
}
//-----------------------------------------------------------------------------
uint HexahedronCell::findEdge(uint i, Cell const& cell) const
{
  // Ordering convention for edges (order of non-incident vertices)

  // Get vertices and edges
  uint const* v = cell.entities(0);
  dolfin_assert(v);
  uint const* e = cell.entities(1);
  dolfin_assert(e);

  // Look for edge satisfying ordering convention
  MeshTopology const& topology = cell.mesh().topology();
  for (uint j = 0; j < 12; ++j)
  {
    uint const * ev = topology(1, 0)(e[j]);
    dolfin_assert(ev);
    // Check incident pairs instead of non-incident quadruples
    uint const v0 = v[EIV[i][0]];
    uint const v1 = v[EIV[i][1]];
    if ((ev[0] == v0 && ev[1] == v1) || (ev[0] == v1 && ev[1] == v0))
    {
      return j;
    }
  }

  // We should not reach this
  error("Unable to find edge with index %d in hexahedron.", cell.index());

  return 0;
}
//-----------------------------------------------------------------------------
uint HexahedronCell::findFace(uint i, Cell const& cell) const
{
  // Ordering convention for faces (order of non-incident vertices)

  // Get vertices and edges
  uint const* v = cell.entities(0);
  dolfin_assert(v);
  uint const* f = cell.entities(2);
  dolfin_assert(f);

  // Look for edge satisfying ordering convention
  MeshTopology const& topology = cell.mesh().topology();
  uint const v0 = v[FNV[i][0]];
  uint const v1 = v[FNV[i][1]];
  uint const v2 = v[FNV[i][2]];
  uint const v3 = v[FNV[i][3]];
  for (uint j = 0; j < 6; ++j)
  {
    uint const * fv = topology(2, 0)(f[j]);
    dolfin_assert(fv);
    if (fv[0] != v0 && fv[0] != v1 && fv[0] != v2 && fv[0] != v3 &&
        fv[1] != v0 && fv[1] != v1 && fv[1] != v2 && fv[1] != v3 &&
        fv[2] != v0 && fv[2] != v1 && fv[2] != v2 && fv[2] != v3 &&
        fv[3] != v0 && fv[3] != v1 && fv[3] != v2 && fv[3] != v3)
    {
      return j;
    }
  }

  // We should not reach this
  error("Unable to find face with index %d in hexahedron.", cell.index());

  return 0;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
