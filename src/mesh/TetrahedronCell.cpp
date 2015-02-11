// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman 2006.
// Modified by Garth N. Wells 2006.
// Modified by Kristian Oelgaard 2006.
//
// First added:  2006-06-05
// Last changed: 2008-06-20

#include <algorithm>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/TetrahedronCell.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/GeometricPredicates.h>
#include <dolfin/parameter/parameters.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
TetrahedronCell::TetrahedronCell() :
    CellType(tetrahedron, triangle)
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
  switch (dim)
    {
    case 0:
      return 4;  // vertices
    case 1:
      return 6;  // edges
    case 2:
      return 4;  // faces
    case 3:
      return 1;  // cells
    default:
      error("Illegal topological dimension %d for tetrahedron.", dim);
      break;
    }

  return 0;
}
//-----------------------------------------------------------------------------
uint TetrahedronCell::numVertices(uint dim) const
{
  switch (dim)
    {
    case 0:
      return 1;  // vertices
    case 1:
      return 2;  // edges
    case 2:
      return 3;  // faces
    case 3:
      return 4;  // cells
    default:
      error("Illegal topological dimension %d for tetrahedron.", dim);
      break;
    }

  return 0;
}
//-----------------------------------------------------------------------------
uint TetrahedronCell::orientation(Cell const& cell) const
{
  // Check that we get a tetrahedron
  dolfin_assert(cell.dim() == 3);
  dolfin_assert(cell.numEntities(0) == 4);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = cell.mesh().geometry();
  uint const * vertices = cell.entities(0);
  real const * v0 = geometry.x(vertices[0]);
  real const * v1 = geometry.x(vertices[1]);
  real const * v2 = geometry.x(vertices[2]);
  real const * v3 = geometry.x(vertices[3]);

  real a =
      + ((v1[1] - v0[1]) * (v2[2] - v0[2]) - (v1[2] - v0[2]) * (v2[1] - v0[1]))
          * (v3[0] - v0[0])
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
      error("Don't know how to create entities of topological dimension %d.",
            dim);
      break;
    }
}
//-----------------------------------------------------------------------------
void TetrahedronCell::orderEntities(Cell& cell) const
{
  // Sort i - j for i > j: 1 - 0, 2 - 0, 2 - 1, 3 - 0, 3 - 1, 3 - 2

  // Check that we get a tetrahedron
  dolfin_assert(cell.dim() == 3);
  dolfin_assert(cell.numEntities(0) == 4);

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

  // Sort local edges on cell after non-incident vertex tuble, connectivity 3-1
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
void TetrahedronCell::refineCell(Cell& cell, MeshEditor& editor,
                                 uint& current_cell) const
{
  // Check that we get a tetrahedron
  dolfin_assert(cell.dim() == 3);
  dolfin_assert(cell.numEntities(0) == 4);

  // Get vertices and edges
  uint const* v = cell.entities(0);
  uint const* e = cell.entities(1);
  dolfin_assert(v);
  dolfin_assert(e);

  // Get offset for new vertex indices
  uint const offset = cell.mesh().numVertices();

  // Compute indices for the ten new vertices
  uint const v0 = v[0];
  uint const v1 = v[1];
  uint const v2 = v[2];
  uint const v3 = v[3];
  uint const e0 = offset + e[findEdge(0, cell)];
  uint const e1 = offset + e[findEdge(1, cell)];
  uint const e2 = offset + e[findEdge(2, cell)];
  uint const e3 = offset + e[findEdge(3, cell)];
  uint const e4 = offset + e[findEdge(4, cell)];
  uint const e5 = offset + e[findEdge(5, cell)];

  // Regular refinement: 8 new cells
  uint connectivity[32] = { v0, e3, e4, e5, v1, e1, e2, e5, v2, e0, e2, e4, v3,
                            e0, e1, e3, e0, e1, e2, e5, e0, e1, e3, e5, e0, e2,
                            e4, e5, e0, e3, e4, e5 };
  for (uint i = 0; i < 8; ++i)
  {
    editor.addCell(current_cell++, &connectivity[i * 4]);
  }
}
//-----------------------------------------------------------------------------
void TetrahedronCell::refineCellIrregular(Cell& cell, MeshEditor& editor,
                                          uint& current_cell,
                                          uint refinement_rule,
                                          uint* marked_edges) const
{
  error("Not implemented yet.");

  /*
   // Get vertices and edges
   uint const* v = cell.entities(0);
   uint const* e = cell.entities(1);
   dolfin_assert(v);
   dolfin_assert(e);

   // Get offset for new vertex indices
   uint const offset = cell.mesh().numVertices();

   // Compute indices for the ten new vertices
   uint const v0 = v[0];
   uint const v1 = v[1];
   uint const v2 = v[2];
   uint const v3 = v[3];
   uint const e0 = offset + e[0];
   uint const e1 = offset + e[1];
   uint const e2 = offset + e[2];
   uint const e3 = offset + e[3];
   uint const e4 = offset + e[4];
   uint const e5 = offset + e[5];

   // Refine according to refinement rule
   // The rules are numbered according to the paper:
   // J. Bey, "Tetrahedral Grid Refinement", 1995.
   switch ( refinement_rule )
   {
   case 1:
   // Rule 1: 4 new cells
   editor.addCell(current_cell++, v0, e1, e3, e2);
   editor.addCell(current_cell++, v1, e2, e4, e0);
   editor.addCell(current_cell++, v2, e0, e5, e1);
   editor.addCell(current_cell++, v3, e5, e4, e3);
   break;
   case 2:
   // Rule 2: 2 new cells
   editor.addCell(current_cell++, v0, e1, e3, e2);
   editor.addCell(current_cell++, v1, e2, e4, e0);
   break;
   case 3:
   // Rule 3: 3 new cells
   editor.addCell(current_cell++, v0, e1, e3, e2);
   editor.addCell(current_cell++, v1, e2, e4, e0);
   editor.addCell(current_cell++, v2, e0, e5, e1);
   break;
   case 4:
   // Rule 4: 4 new cells
   editor.addCell(current_cell++, v0, e1, e3, e2);
   editor.addCell(current_cell++, v1, e2, e4, e0);
   editor.addCell(current_cell++, v2, e0, e5, e1);
   editor.addCell(current_cell++, v3, e5, e4, e3);
   break;
   default:
   error("Illegal rule for irregular refinement of tetrahedron.");
   }
   */
}
//-----------------------------------------------------------------------------
real TetrahedronCell::volume(MeshEntity const& tetrahedron) const
{
  // Check that we get a tetrahedron
  dolfin_assert(tetrahedron.dim() == 3);
  dolfin_assert(tetrahedron.numEntities(0) == 4);

  // Get mesh geometry
  MeshGeometry const& geometry = tetrahedron.mesh().geometry();

  // Only know how to compute the volume when embedded in R^3
  dolfin_assert(geometry.dim() > 2);

  // Get the coordinates of the four vertices
  uint const* vertices = tetrahedron.entities(0);
  real const* x0 = geometry.x(vertices[0]);
  real const* x1 = geometry.x(vertices[1]);
  real const* x2 = geometry.x(vertices[2]);
  real const* x3 = geometry.x(vertices[3]);

  // Formula for volume from http://mathworld.wolfram.com
  real V = (
      + x0[0]
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
real TetrahedronCell::diameter(MeshEntity const& tetrahedron) const
{
  // Check that we get a tetrahedron
  dolfin_assert(tetrahedron.dim() == 3);
  dolfin_assert(tetrahedron.numEntities(0) == 4);

  // Get mesh geometry
  MeshGeometry const& geometry = tetrahedron.mesh().geometry();

  // Only know how to compute the volume when embedded in R^3
  dolfin_assert(geometry.dim() > 2);

  // Get the coordinates of the four vertices
  uint const* vertices = tetrahedron.entities(0);
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
real TetrahedronCell::circumradius(MeshEntity const& tetrahedron) const
{
  // Check that we get a tetrahedron
  dolfin_assert(tetrahedron.dim() == 3);
  dolfin_assert(tetrahedron.numEntities(0) == 4);

  // Get mesh geometry
  MeshGeometry const& geometry = tetrahedron.mesh().geometry();

  // Only know how to compute the volume when embedded in R^3
  dolfin_assert(geometry.dim() > 2);

  // Get the coordinates of the four vertices
  uint const* vertices = tetrahedron.entities(0);
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
  real area = sqrt(s * (s - la) * (s - lb) * (s - lc));

  // Formula for volume from http://mathworld.wolfram.com
  real V = (
        + x0[0]
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
real TetrahedronCell::normal(Cell const& cell, uint facet, uint i) const
{
  return normal(cell, facet)[i];
}
//-----------------------------------------------------------------------------
Point TetrahedronCell::normal(Cell const& cell, uint facet) const
{
  // This is a trick to be allowed to initialize a facet from the cell
  Cell& c = const_cast<Cell&>(cell);

  // Create facet from the mesh and local facet number
  Facet f(c.mesh(), c.entities(2)[facet]);

  // Get global index of opposite vertex
  uint const v0 = cell.entities(0)[facet];

  // Get global index of vertices on the facet
  uint v1 = f.entities(0)[0];
  uint v2 = f.entities(0)[1];
  uint v3 = f.entities(0)[2];

  // Get mesh geometry
  MeshGeometry const& geometry = cell.mesh().geometry();

  // Get the coordinates of the four vertices
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
  dolfin_assert(cell.mesh().topology().dim() == 3);
  dolfin_assert(cell.mesh().geometry().dim() == 3);

  // This is a trick to be allowed to initialize a facet from the cell
  Cell& c = const_cast<Cell&>(cell);

  // Create facet from the mesh and local facet number
  Facet f(c.mesh(), c.entities(2)[facet]);

  // Get mesh geometry
  MeshGeometry const& geometry = cell.mesh().geometry();

  // Get the coordinates of the three vertices
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
  return 0.5 * sqrt(v0 * v0 + v1 * v1 + v2 * v2);
}
//-----------------------------------------------------------------------------
bool TetrahedronCell::intersects(MeshEntity const& tetrahedron,
                                 Point const& p) const
{
  // Adapted from gts_point_is_in_triangle from GTS

  // Get mesh geometry
  MeshGeometry const& geometry = tetrahedron.mesh().geometry();

  // Get global index of vertices of the tetrahedron
  uint v0 = tetrahedron.entities(0)[0];
  uint v1 = tetrahedron.entities(0)[1];
  uint v2 = tetrahedron.entities(0)[2];
  uint v3 = tetrahedron.entities(0)[3];

  // Check orientation
  uint vtmp;
  if (orientation((Cell&) tetrahedron) == 1)
  {
    vtmp = v3;
    v3 = v2;
    v2 = vtmp;
  }

  // Get the coordinates of the four vertices
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
bool TetrahedronCell::intersects(MeshEntity const& interval, Point const& p1,
                                 Point const& p2) const
{
  // FIXME: Not implemented
  error("TetrahedronCell::intersects() not implemented");

  return false;
}
//-----------------------------------------------------------------------------
std::string TetrahedronCell::description() const
{
  std::string s = "tetrahedron (simplex of topological dimension 3)";
  return s;
}
//-----------------------------------------------------------------------------
uint TetrahedronCell::findEdge(uint i, Cell const& cell) const
{
  // Get vertices and edges
  uint const* v = cell.entities(0);
  uint const* e = cell.entities(1);
  dolfin_assert(v);
  dolfin_assert(e);

  // Ordering convention for edges (order of non-incident vertices)
  static uint EV[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 }, {
      2, 3 } };

  // Look for edge satisfying ordering convention
  for (uint j = 0; j < 6; ++j)
  {
    uint const* ev = cell.mesh().topology()(1, 0)(e[j]);
    dolfin_assert(ev);
    uint const v0 = v[EV[i][0]];
    uint const v1 = v[EV[i][1]];
    if (ev[0] != v0 && ev[0] != v1 && ev[1] != v0 && ev[1] != v1) return j;
  }

  // We should not reach this
  error("Unable to find edge in tetrahedron cell with index %d.", cell.index());

  return 0;
}
//-----------------------------------------------------------------------------

}
