// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
//

#include <dolfin/mesh/HexahedronCell.h>

#include <dolfin/common/maybe_unused.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshEditor.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------

// UFC: Number of Entities
uint const HexahedronCell::NE[4][4] =
{ { 1, 0, 0, 0 }, { 2, 1, 0, 0 }, { 4, 4, 1, 0 }, { 8, 12, 6, 1 } };

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
    CellType("hexahedron", CellType::hexahedron, CellType::quadrilateral)
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
  return NE[3][dim];
}
//-----------------------------------------------------------------------------
uint HexahedronCell::num_entities(uint d0, uint d1) const
{
  dolfin_assert(d0 <= TD);
  dolfin_assert(d1 <= TD);
  return NE[d0][d1];
}
//-----------------------------------------------------------------------------
uint HexahedronCell::num_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NE[dim][0];
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
      e[1][2] = v[7];
      e[1][3] = v[6];
      e[2][0] = v[1];
      e[2][1] = v[2];
      e[2][2] = v[6];
      e[2][3] = v[5];
      e[3][0] = v[0];
      e[3][1] = v[4];
      e[3][2] = v[7];
      e[3][3] = v[3];
      e[4][0] = v[0];
      e[4][1] = v[1];
      e[4][2] = v[5];
      e[4][3] = v[4];
      e[5][0] = v[0];
      e[5][1] = v[3];
      e[5][2] = v[2];
      e[5][3] = v[1];
      break;
    default:
      error("Invalid topological dimension for creation of entities: %d.", dim);
      break;
    }
}
//-----------------------------------------------------------------------------
void HexahedronCell::order_entities(MeshTopology& topology, uint i) const
{
  // Sort i - j for i > j: 1 - 0, 2 - 0, 2 - 1, 3 - 0, 3 - 1, 3 - 2
  dolfin_assert(topology.type(i).cellType() == this->cell_type);

  // Sort local vertices on edges in ascending order, connectivity 1 - 0
  if (topology.connectivity(1, 0))
  {
    dolfin_assert(topology.connectivity(3, 1));

    // Get edges
    uint* cell_edges = topology(3, 1)(i);

    // Sort vertices on each edge
    for (uint i = 0; i < 12; ++i)
    {
      uint* edge_vertices = topology(1, 0)(cell_edges[i]);
      std::sort(edge_vertices, edge_vertices + 2);
    }
  }

  // Sort local vertices on facets in ascending order, connectivity 2 - 0
  if (topology.connectivity(2, 0))
  {
    dolfin_assert(topology.connectivity(3, 2));

    // Get facets
    uint* cell_facets = topology(3, 2)(i);

    // Sort vertices on each facet
    for (uint i = 0; i < 6; ++i)
    {
      uint* facet_vertices = topology(2, 0)(cell_facets[i]);
      std::sort(facet_vertices, facet_vertices + 4);
    }
  }

  // Sort local edges on local facets after non-incident vertex, connectivity 2 - 1
  if (topology.connectivity(2, 1))
  {
    dolfin_assert(topology.connectivity(3, 2));
    dolfin_assert(topology.connectivity(2, 0));
    dolfin_assert(topology.connectivity(1, 0));

    // Get facet numbers
    uint* cell_facets = topology(2, 1)(i);

    // Loop over facets on cell
    for (uint i = 0; i < 6; ++i)
    {
      // For each facet number get the global vertex numbers
      // Facet
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

  // Vertices are supposed to have been ordered appropriately

  // Sort local edges on cell after non-incident vertex tuple, connectivity 3-1
  if (topology.connectivity(3, 1))
  {
    dolfin_assert(topology.connectivity(1, 0));

    // Get cell vertices and edge numbers
    uint* cell_vertices = topology(3, 0)(i);
    uint* cell_edges = topology(3, 1)(i);

    // Loop two vertices on cell as a lexicographical tuple
    uint m = 0;
    for (uint t = 0; t < 12; ++t)
    {
      uint v0 = ENV[t][0];
      uint v1 = ENV[t][1];
      // Loop edge numbers
      for (uint k = m; k < 12; ++k)
      {
        // Get local vertices on edge
        uint* edge_vertices = topology(1, 0)(cell_edges[k]);

        // Check if the ith and jth vertex of the cell are non-incident on edge k
#if __SUNPRO_CC
        int n1 = 0;
        int n2 = 0;
        std::count(edge_vertices, edge_vertices + 2, cell_vertices[v0], n1);
        std::count(edge_vertices, edge_vertices + 2, cell_vertices[v1], n2);
        if (!n1 && !n2 )
#else
        if (!std::count(edge_vertices, edge_vertices + 2, cell_vertices[v0])
         && !std::count(edge_vertices, edge_vertices + 2, cell_vertices[v1]))
#endif
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

  // Sort local facets on cell after non-incident vertex, connectivity 3 - 2
  if (topology.connectivity(3, 2))
  {
    dolfin_assert(topology.connectivity(2, 0));

    // Get cell vertices and facet numbers
    uint* cell_vertices = topology(3, 0)(i);
    uint* cell_facets = topology(3, 2)(i);

    //
    uint m = 0;
    for (uint t = 0; t < 6; ++t)
    {
      uint v0 = FNV[t][0];
      uint v1 = FNV[t][1];
      uint v2 = FNV[t][2];
      uint v3 = FNV[t][3];

      // Loop facets on cell
      for (uint k = m; k < 6; ++k)
      {
        uint* facet_vertices = topology(2, 0)(cell_facets[k]);

        // Check if the ith vertex of the cell is non-incident on facet j
#if __SUNPRO_CC
        int n1 = 0;
        int n2 = 0;
        int n3 = 0;
        int n4 = 0;
        std::count(facet_vertices, facet_vertices + 4, cell_vertices[v0], n1);
        std::count(facet_vertices, facet_vertices + 4, cell_vertices[v1], n2);
        std::count(facet_vertices, facet_vertices + 4, cell_vertices[v2], n3);
        std::count(facet_vertices, facet_vertices + 4, cell_vertices[v3], n4);
        if (!n1 && !n2 && !n3 && !n4)
#else
        if (!std::count(facet_vertices, facet_vertices + 4, cell_vertices[v0])
         && !std::count(facet_vertices, facet_vertices + 4, cell_vertices[v1])
         && !std::count(facet_vertices, facet_vertices + 4, cell_vertices[v2])
         && !std::count(facet_vertices, facet_vertices + 4, cell_vertices[v3]))
#endif
        {
          // Swap facet numbers
          uint tmp = cell_facets[m];
          cell_facets[m] = cell_facets[k];
          cell_facets[k] = tmp;
          ++m;
          break;
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
void HexahedronCell::order_facet(uint vertices[], Facet& facet) const
{
  static uint e[4] = { 0 };
  Mesh& mesh = facet.mesh();
  Cell const cell(mesh, facet.entities(TD)[0]);

  // Find facet, it is outward oriented by definition.
  //NOTE: Keep consistency with create_entities.
  uint facet_index = cell.index(facet);
  uint const * v = cell.entities(0);
  switch(facet_index)
  {
    case 0:
      e[0] = v[4];
      e[1] = v[5];
      e[2] = v[6];
      e[3] = v[7];
      break;
    case 1:
      e[0] = v[2];
      e[1] = v[3];
      e[2] = v[7];
      e[3] = v[6];
      break;
    case 2:
      e[0] = v[1];
      e[1] = v[2];
      e[2] = v[6];
      e[3] = v[5];
      break;
    case 3:
      e[0] = v[0];
      e[1] = v[4];
      e[2] = v[7];
      e[3] = v[3];
      break;
    case 4:
      e[0] = v[0];
      e[1] = v[1];
      e[2] = v[5];
      e[3] = v[4];
      break;
    case 5:
      e[0] = v[0];
      e[1] = v[3];
      e[2] = v[2];
      e[3] = v[1];
      break;
    default:
      error("HexahedronCell : invalid local facet index %u", facet_index);
      break;
  }

  // Reorder facet vertices
  for (uint i = 0; i < 4; ++i)
  {
    for (uint j = 0; j < 4; ++j)
    {
      if(e[i] == facet.entities(0)[j])
      {
        e[i] = vertices[j];
        break;
      }
    }
  }
  std::memcpy(&vertices[0], &e[0], 4 * sizeof(uint));
}
//-----------------------------------------------------------------------------
bool HexahedronCell::connectivity_needs_ordering(uint d0, uint d1) const
{
  dolfin_assert(d0 <= TD && d1 <= TD);
  // Do not order cell - vertices connectivities
  return (d0 > 0 && d0 > d1) && !(d0 == TD && d1 == 0);
}
//-----------------------------------------------------------------------------
void HexahedronCell::initialize_connectivities(Mesh& mesh) const
{
  mesh.init(1, 0);
  mesh.init(2, 0);
  mesh.init(2, 1);
  mesh.init(3, 0);
  mesh.init(3, 1);
  mesh.init(3, 2);
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
  uint const eoffset = cell.mesh().size(0);
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
  uint const foffset = eoffset + cell.mesh().size(1);
  uint const f00 = foffset + f[findFace(0, cell)];
  uint const f01 = foffset + f[findFace(1, cell)];
  uint const f02 = foffset + f[findFace(2, cell)];
  uint const f03 = foffset + f[findFace(3, cell)];
  uint const f04 = foffset + f[findFace(4, cell)];
  uint const f05 = foffset + f[findFace(5, cell)];
  uint const coffset = foffset + cell.mesh().size(2);
  uint const c00 = coffset + cell.index();

  // Add the eight new cells
  uint const cv0[8] = { v00, e11, f05, e10, e09, f04, c00, f03 };
  editor.add_cell(current_cell++, &cv0[0]);
  uint const cv1[8] = { e11, v01, e08, f05, f04, e07, f02, c00 };
  editor.add_cell(current_cell++, &cv1[0]);
  uint const cv2[8] = { f05, e08, v02, e06, c00, f02, e05, f01 };
  editor.add_cell(current_cell++, &cv2[0]);
  uint const cv3[8] = { e10, f05, e06, v03, f03, c00, f01, e04 };
  editor.add_cell(current_cell++, &cv3[0]);
  uint const cv4[8] = { e09, f04, c00, f03, v04, e03, f00, e02 };
  editor.add_cell(current_cell++, &cv4[0]);
  uint const cv5[8] = { f04, e07, f02, c00, e03, v05, e01, f00 };
  editor.add_cell(current_cell++, &cv5[0]);
  uint const cv6[8] = { c00, f02, e05, f01, f00, e01, v06, e00 };
  editor.add_cell(current_cell++, &cv6[0]);
  uint const cv7[8] = { f03, c00, f01, e04, e02, f00, e00, v07 };
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
  MAYBE_UNUSED(dim);
  return 1;
}
//-----------------------------------------------------------------------------
real HexahedronCell::volume(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[3][0]);

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
  dolfin_assert(entity.num_entities(0) == NE[3][0]);

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
  //FIXME: No better idea for now
  return 0.5*std::sqrt(3.0)*this->diameter(entity);
}
//-----------------------------------------------------------------------------
real HexahedronCell::inradius(MeshEntity const& entity) const
{
  //FIXME: No better idea for now
  return 0.5*this->diameter(entity);
}
//-----------------------------------------------------------------------------
void HexahedronCell::midpoint(MeshEntity const& entity, real * p) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[3][0]);

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
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    p[i] = 0.125
        * (x0[i] + x1[i] + x2[i] + x3[i] + x4[i] + x5[i] + x6[i] + x7[i]);
  }
}
//-----------------------------------------------------------------------------
void HexahedronCell::normal(Cell const& cell, uint facet, real * n) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Create facet from the mesh and local facet number
  Cell& c = const_cast<Cell&>(cell);
  Facet f(c.mesh(), c.entities(1)[facet]);
  MeshGeometry const& geometry = cell.mesh().geometry();

  // first non-incident vertex
  real const * vf = geometry.x(c.entities(0)[FNV[facet][0]]);

  // vertices on the facet
  uint const * vertices = f.entities(0);
  real const * v0 = geometry.x(vertices[FIV[facet][0]]);
  real const * v1 = geometry.x(vertices[FIV[facet][1]]);
  real const * v3 = geometry.x(vertices[FIV[facet][3]]);

  // Vector normal to facet
  n[0] = + (v1[1] - v0[1]) * (v3[2] - v0[2])
         - (v1[2] - v0[2]) * (v3[1] - v0[1]);
  n[1] = + (v1[2] - v0[2]) * (v3[0] - v0[0])
         - (v1[0] - v0[0]) * (v3[2] - v0[2]);
  n[2] = + (v1[0] - v0[0]) * (v3[1] - v0[1])
         - (v1[1] - v0[1]) * (v3[0] - v0[0]);
  real const nn = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  n[0] /= nn;
  n[1] /= nn;
  n[2] /= nn;

  // Flip direction of normal so it points outward
  if (+ n[0]*(vf[0]-v0[0]) + n[1]*(vf[1]-v0[1]) + n[2]*(vf[2] - v0[2]) > 0.0)
  {
    n[0] *= -1.0;
    n[1] *= -1.0;
    n[2] *= -1.0;
  }
}
//-----------------------------------------------------------------------------
real HexahedronCell::facet_area(Cell const& cell, uint facet) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Create facet from the mesh and local facet number
  Cell& c = const_cast<Cell&>(cell);
  Facet f(c.mesh(), c.entities(1)[facet]);

  // Get the coordinates of the two vertices
  MeshGeometry const& geometry = cell.mesh().geometry();
  real const * p = geometry.x(f.entities(0)[FIV[facet][0]]);
  real const * q = geometry.x(f.entities(0)[FIV[facet][1]]);
  real const * r = geometry.x(f.entities(0)[FIV[facet][2]]);
  real const * s = geometry.x(f.entities(0)[FIV[facet][3]]);

  // Compute the area with diagonal formula
  real c0 = (r[1] - p[1]) * (s[2] - q[2]) - (r[2] - p[2]) * (s[1] - q[1]);
  real c1 = (r[2] - p[2]) * (s[0] - q[0]) - (r[0] - p[0]) * (s[2] - q[2]);
  real c2 = (r[0] - p[0]) * (s[1] - q[1]) - (r[1] - p[1]) * (s[0] - q[0]);
  return 0.5 * std::sqrt(c0 * c0 + c1 * c1 + c2 * c2);
}
//-----------------------------------------------------------------------------
bool HexahedronCell::intersects(MeshEntity const& e, Point const&) const
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.num_entities(0) == NE[3][0]);
  MAYBE_UNUSED(e)

  // Get the coordinates of the vertices
  /*
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
   */

  error("Collision of hexahedron with point not implemented.");

  return true;
}
//-----------------------------------------------------------------------------
bool HexahedronCell::intersects(MeshEntity const& e, Point const&,
                                Point const&) const
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.num_entities(0) == NE[3][0]);
  MAYBE_UNUSED(e)

  // Get the coordinates of the vertices
  /*
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
   */

  error("Collision of hexahedron with segment not implemented.");

  return true;
}
//-----------------------------------------------------------------------------
void HexahedronCell::create_reference_cell(Mesh& mesh) const
{
  MeshEditor me(mesh, CellType::hexahedron, 3, DOLFIN_COMM_SELF);
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
}
//-----------------------------------------------------------------------------
real const * HexahedronCell::reference_vertex(uint i) const
{
  return &VC[i][0];
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
bool HexahedronCell::check(Cell& cell) const
{
  bool ret = CellType::check(cell);

  //
  MeshTopology const& topology = cell.mesh().topology();
  uint const* v = cell.entities(0);
  dolfin_assert(v);

  // Check edge -> incident vertices mapping
  if (cell.mesh().topology().connectivity(1, 0))
  {
    uint const* e = cell.entities(1);
    dolfin_assert(e);
    for (uint i = 0; i < 12; ++i)
    {
      uint const * ev = topology(1, 0)(e[i]);
      dolfin_assert(ev);
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

  // Check face -> incident vertices mapping
  if (cell.mesh().topology().connectivity(2, 0))
  {
    uint const* f = cell.entities(2);
    dolfin_assert(f);
    for (uint i = 0; i < 6; ++i)
    {
      uint const * fv = topology(2, 0)(f[i]);
      dolfin_assert(fv);
      for (uint j = 0; j < 4; ++j)
      {
        if (fv[j] != v[FIV[i][0]] && fv[j] != v[FIV[i][1]] &&
            fv[j] != v[FIV[i][2]] && fv[j] != v[FIV[i][3]])
        {
          ret = false;
          warning("CellType : invalid face -> incident vertices mapping\n"
                "f[%u] = (v%u, v%u, v%u, v%u)", i, fv[0], fv[1], fv[2], fv[3]);
        }
      }
    }
  }

  return ret;
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
