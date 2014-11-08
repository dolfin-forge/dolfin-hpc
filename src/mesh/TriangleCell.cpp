// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2006.
// Modified by Kristian Oelgaard, 2006-2007.
// Modified by Dag Lindbo, 2008
// Modified by Aurelien Larcher, 2014
//
// First added:  2006-06-05
// Last changed: 2014-11-07

#include <dolfin/mesh/TriangleCell.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/GeometricPredicates.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/parameter/parameters.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------
TriangleCell::TriangleCell() :
    CellType(triangle, interval)
{
}

//-----------------------------------------------------------------------------
uint TriangleCell::dim() const
{
  return 2;
}
//-----------------------------------------------------------------------------
uint TriangleCell::numEntities(uint dim) const
{
  switch (dim)
    {
    case 0:
      return 3;  // vertices
    case 1:
      return 3;  // edges
    case 2:
      return 1;  // cells
    default:
      error("Illegal topological dimension %d for triangle.", dim);
      break;
    }

  return 0;
}
//-----------------------------------------------------------------------------
uint TriangleCell::numVertices(uint dim) const
{
  switch (dim)
    {
    case 0:
      return 1;  // vertices
    case 1:
      return 2;  // edges
    case 2:
      return 3;  // cells
    default:
      error("Illegal topological dimension %d for triangle.", dim);
      break;
    }

  return 0;
}
//-----------------------------------------------------------------------------
uint TriangleCell::orientation(Cell const& cell) const
{
  // Check that we get a triangle
  dolfin_assert(cell.dim() == 2);
  dolfin_assert(cell.numEntities(0) == 3);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = cell.mesh().geometry();
  uint const * vertices = cell.entities(0);
  real const * v0 = geometry.x(vertices[0]);
  real const * v1 = geometry.x(vertices[1]);
  real const * v2 = geometry.x(vertices[2]);
  return (
      ((v1[0] - v0[0]) * (v2[1] - v0[1]) - (v1[1] - v0[1]) * (v2[0] - v0[0]))
          < 0.0 ? 1 : 0);
}
//-----------------------------------------------------------------------------
void TriangleCell::createEntities(uint** e, uint dim, uint const* v) const
{
  // We only need to know how to create edges
  if (dim != 1)
  {
    error("Don't know how to create entities of topological dimension %d.",
          dim);
  }

  // Create the three edges
  e[0][0] = v[1];
  e[0][1] = v[2];
  e[1][0] = v[0];
  e[1][1] = v[2];
  e[2][0] = v[0];
  e[2][1] = v[1];
}
//-----------------------------------------------------------------------------
void TriangleCell::orderEntities(Cell& cell) const
{
  // Sort i - j for i > j: 1 - 0, 2 - 0, 2 - 1

  // Get mesh topology
  MeshTopology& topology = cell.mesh().topology();

  // Sort local vertices on edges in ascending order, connectivity 1 - 0
  if (topology(1, 0).size() > 0)
  {
    dolfin_assert(topology(2, 1).size() > 0);

    // Get edges
    uint* cell_edges = cell.entities(1);

    // Sort vertices on each edge
    for (uint i = 0; i < 3; ++i)
    {
      uint* edge_vertices = topology(1, 0)(cell_edges[i]);
      std::sort(edge_vertices, edge_vertices + 2);
    }
  }

  // Sort local vertices on cell in ascending order, connectivity 2 - 0
  if (topology(2, 0).size() > 0)
  {
    uint* cell_vertices = cell.entities(0);
    std::sort(cell_vertices, cell_vertices + 3);
  }

  // Sort local edges on cell after non-incident vertex, connectivity 2 - 1
  if (topology(2, 1).size() > 0)
  {
    dolfin_assert(topology(2, 1).size() > 0);

    // Get cell vertices and edges
    uint* cell_vertices = cell.entities(0);
    uint* cell_edges = cell.entities(1);

    // Loop over vertices on cell
    for (uint i = 0; i < 3; ++i)
    {
      // Loop over edges on cell
      for (uint j = i; j < 3; ++j)
      {
        uint* edge_vertices = topology(1, 0)(cell_edges[j]);

        // Check if the ith vertex of the cell is non-incident with edge j
#if __SUNPRO_CC
        int n1 = 0;
        std::count(edge_vertices, edge_vertices + 2, cell_vertices[i], n1);
        if ( n1 == 0)
#else
        if (std::count(edge_vertices, edge_vertices + 2, cell_vertices[i]) == 0)
#endif
        {
          // Swap edge numbers
          uint tmp = cell_edges[i];
          cell_edges[i] = cell_edges[j];
          cell_edges[j] = tmp;
          break;
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
void TriangleCell::refineCell(Cell& cell, MeshEditor& editor,
                              uint& current_cell) const
{
  // Check that we get a triangle
  dolfin_assert(cell.dim() == 2);
  dolfin_assert(cell.numEntities(0) == 3);

  // Get vertices and edges
  uint const * v = cell.entities(0);
  uint const * e = cell.entities(1);
  dolfin_assert(v);
  dolfin_assert(e);

  // Get offset for new vertex indices
  uint const offset = cell.mesh().numVertices();

  // Compute indices for the six new vertices
  uint const v0 = v[0];
  uint const v1 = v[1];
  uint const v2 = v[2];
  uint const e0 = offset + e[findEdge(0, cell)];
  uint const e1 = offset + e[findEdge(1, cell)];
  uint const e2 = offset + e[findEdge(2, cell)];

  // Add the four new cells
  uint connectivity[12] = { v0, e2, e1, v1, e0, e2, v2, e1, e0, e0, e1, e2 };
  for (uint i = 0; i < 4; ++i)
  {
    editor.addCell(current_cell++, &connectivity[i * 3]);
  }
}
//-----------------------------------------------------------------------------
real TriangleCell::volume(MeshEntity const& triangle) const
{
  // Check that we get a triangle
  dolfin_assert(triangle.dim() == 2);
  dolfin_assert(triangle.numEntities(0) == 3);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = triangle.mesh().geometry();
  uint const * vertices = triangle.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);

  switch (geometry.dim())
    {
    case 2:
      // Compute area of triangle embedded in R^2
      // Formula for volume from http://mathworld.wolfram.com
      return 0.5
          * std::abs(
              (x0[0] * x1[1] + x0[1] * x2[0] + x1[0] * x2[1])
                  - (x2[0] * x1[1] + x2[1] * x0[0] + x1[0] * x0[1]));
      break;
    case 3:
      // Compute area of triangle embedded in R^3
      // Formula for volume from http://mathworld.wolfram.com
      return 0.5 * std::sqrt(
      + std::pow( (x0[1] * x1[2] + x0[2] * x2[1] + x1[1] * x2[2])
                    - (x2[1] * x1[2] + x2[2] * x0[1] + x1[1] * x0[2]), 2)
      + std::pow( (x0[2] * x1[0] + x0[0] * x2[2] + x1[2] * x2[0])
                    - (x2[2] * x1[0] + x2[0] * x0[2] + x1[2] * x0[0]), 2)
      + std::pow( (x0[0] * x1[1] + x0[1] * x2[0] + x1[0] * x2[1])
                    - (x2[0] * x1[1] + x2[1] * x0[0] + x1[0] * x0[1]), 2)
      );
      break;
    default:
      error("Implemented measure of a triangle only embedded in R^2 or R^3.");
      break;
    }
  return 0.0;
}
//-----------------------------------------------------------------------------
real TriangleCell::diameter(MeshEntity const& triangle) const
{
  // Check that we get a triangle
  dolfin_assert(triangle.dim() == 2);
  dolfin_assert(triangle.numEntities(0) == 3);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = triangle.mesh().geometry();
  uint const * vertices = triangle.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real e0 = 0.0;
  real e1 = 0.0;
  real e2 = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    e0 += (x1[i] - x0[i]) * (x1[i] - x0[i]);
    e1 += (x2[i] - x1[i]) * (x2[i] - x1[i]);
    e2 += (x0[i] - x2[i]) * (x0[i] - x2[i]);
  }
  return std::sqrt(std::max(std::max(e0, e1), e2));
}
//-----------------------------------------------------------------------------
real TriangleCell::circumradius(MeshEntity const& triangle) const
{
  // Check that we get a triangle
  dolfin_assert(triangle.dim() == 2);
  dolfin_assert(triangle.numEntities(0) == 3);

  // Get the coordinates of the three vertices
  MeshGeometry const& geometry = triangle.mesh().geometry();
  uint const * vertices = triangle.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  real const * x2 = geometry.x(vertices[2]);
  real e0 = 0.0;
  real e1 = 0.0;
  real e2 = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    e0 += (x1[i] - x0[i]) * (x1[i] - x0[i]);
    e1 += (x2[i] - x1[i]) * (x2[i] - x1[i]);
    e2 += (x0[i] - x2[i]) * (x0[i] - x2[i]);
  }
  e0 = std::sqrt(e0);
  e1 = std::sqrt(e1);
  e2 = std::sqrt(e2);

  // Formula for circumradius from http://mathworld.wolfram.com
  // Using Heron's formula for the volume instead of calling volume()
  real const s = 0.5 * (e0 + e1 + e2);
  return 0.25 * e0 * e1 * e2 / std::sqrt(s * (s - e0) * (s - e1) * (s - e2));
}
//-----------------------------------------------------------------------------
real TriangleCell::normal(Cell const& cell, uint facet, uint i) const
{
  return normal(cell, facet)[i];
}
//-----------------------------------------------------------------------------
Point TriangleCell::normal(Cell const& cell, uint facet) const
{
  // Check that we get a triangle
  dolfin_assert(cell.dim() == 2);
  dolfin_assert(cell.numEntities(0) == 3);

  // This is a trick to be allowed to initialize a facet from the cell
  Cell& c = const_cast<Cell&>(cell);

  // Get geometry
  MeshGeometry const& geometry = c.mesh().geometry();

  // The normal vector is currently only defined for a triangle in R^2
  if (geometry.dim() != 2)
  {
    error("The normal vector is only defined when the triangle is in R^2");
  }

  // Create facet from the mesh and local facet number
  Facet f(c.mesh(), c.entities(1)[facet]);

  // Get global index of opposite vertex
  uint const v0 = c.entities(0)[facet];

  // Get global index of vertices on the facet
  uint const v1 = f.entities(0)[0];
  uint const v2 = f.entities(0)[1];

  // Get the coordinates of the three vertices
  real const * p0 = geometry.x(v0);
  real const * p1 = geometry.x(v1);
  real const * p2 = geometry.x(v2);

  // Vector normal to facet
  Point n(p2[1] - p1[1], p1[0] - p2[0]);

  // Normalize
  n /= std::sqrt(n[0] * n[0] + n[1] * n[1]);

  // Flip direction of normal so it points outward
  if ((n[0] * (p0[0] - p1[0]) + n[1] * (p0[1] - p1[1])) > 0) n *= -1.0;

  return n;
}
//-----------------------------------------------------------------------------
real TriangleCell::facetArea(Cell const& cell, uint facet) const
{
  // This is a trick to be allowed to initialize a facet from the cell
  Cell& c = const_cast<Cell&>(cell);

  // Create facet from the mesh and local facet number
  Facet f(c.mesh(), c.entities(1)[facet]);

  // Get global index of vertices on the facet
  uint const v0 = f.entities(0)[0];
  uint const v1 = f.entities(0)[1];

  // Get mesh geometry
  MeshGeometry const& geometry = cell.mesh().geometry();

  // Get the coordinates of the two vertices
  real const * p0 = geometry.x(v0);
  real const * p1 = geometry.x(v1);

  // Compute distance between vertices
  real d = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    d += (p0[i] - p1[i]) * (p0[i] - p1[i]);
  }
  return std::sqrt(d);
}
//-----------------------------------------------------------------------------
bool TriangleCell::intersects(MeshEntity const& triangle, Point const& p) const
{
  // Adapted from gts_point_is_in_triangle from GTS

  // Get mesh geometry
  MeshGeometry const& geometry = triangle.mesh().geometry();

  // Get global index of vertices of the triangle
  uint v0 = triangle.entities(0)[0];
  uint v1 = triangle.entities(0)[1];
  uint v2 = triangle.entities(0)[2];

  // Check orientation
  uint vtmp;
  if (orientation((Cell&) triangle) == 1)
  {
    vtmp = v2;
    v2 = v1;
    v1 = vtmp;
  }

  // Get the coordinates of the three vertices
  real const * x0 = geometry.x(v0);
  real const * x1 = geometry.x(v1);
  real const * x2 = geometry.x(v2);

  // Test orientation of p w.r.t. each edge
  real tol = geometry.abs_tolerance(2);
  real d1 = orient2d( x0, x1, &p[0]);
  if (d1 < (-tol)) return false;
  real d2 = orient2d( x1, x2, &p[0]);
  if (d2 < (-tol)) return false;
  real d3 = orient2d( x2, x0, &p[0]);
  if (d3 < (-tol)) return false;

  return true;
}
//-----------------------------------------------------------------------------
bool TriangleCell::intersects(MeshEntity const& tri, Point const& p1,
                              Point const& p2) const
{
  // Adapted from gts_point_is_in_triangle from GTS

  // Get mesh geometry
  MeshGeometry const& geometry = tri.mesh().geometry();

  // Get global index of vertices of the triangle
  uint v0 = tri.entities(0)[0];
  uint v1 = tri.entities(0)[1];
  uint v2 = tri.entities(0)[2];

  // Check orientation
  uint vtmp;
  if (orientation((Cell&) tri) == 1)
  {
    vtmp = v2;
    v2 = v1;
    v1 = vtmp;
  }

  // Get the coordinates of the three vertices
  real const * x0 = geometry.x(v0);
  real const * x1 = geometry.x(v1);
  real const * x2 = geometry.x(v2);

  real d1, d2, d3;

  // Test orientation of each vertex w.r.t. p1-p2
  d1 = orient2d(&p1[0], &p2[0], x0);
  d2 = orient2d(&p1[0], &p2[0], x1);
  d3 = orient2d(&p1[0], &p2[0], x2);

  if (d1 < 0 && d2 < 0 && d3 < 0) return false;
  if (d1 > 0 && d2 > 0 && d3 > 0) return false;

  // Line p1-p2 intersects triangle but both p1 and p2 are
  // on the negative side of x0-x1:
  d1 = orient2d(x0, x1, &p1[0]);
  d2 = orient2d(x0, x1, &p2[0]);

  if (d1 < 0 && d2 < 0) return false;

  // Line p1-p2 intersects triangle but both p1 and p2 are
  // on the negative side of x1-x2:
  d1 = orient2d(x1, x2, &p1[0]);
  d2 = orient2d(x1, x2, &p2[0]);

  if (d1 < 0 && d2 < 0) return false;

  // Line p1-p2 intersects triangle but both p1 and p2 are
  // on the negative side of x2-x0:
  d1 = orient2d(x2, x0, &p1[0]);
  d2 = orient2d(x2, x0, &p2[0]);

  if (d1 < 0 && d2 < 0) return false;

  return true;
}
//-----------------------------------------------------------------------------
std::string TriangleCell::description() const
{
  std::string s = "triangle (simplex of topological dimension 2)";
  return s;
}
//-----------------------------------------------------------------------------
uint TriangleCell::findEdge(uint i, Cell const& cell) const
{
  // Get vertices and edges
  uint const* v = cell.entities(0);
  uint const* e = cell.entities(1);
  dolfin_assert(v);
  dolfin_assert(e);

  // Look for edge satisfying ordering convention
  MeshTopology const& topology = cell.mesh().topology();
  for (uint j = 0; j < 3; ++j)
  {
    uint const * ev = topology(1, 0)(e[j]);
    dolfin_assert(ev);
    if (ev[0] != v[i] && ev[1] != v[i])
    {
      return j;
    }
  }

  // We should not reach this
  error("Unable to find edge in triangle cell with index %d.", cell.index());

  return 0;
}
//-----------------------------------------------------------------------------

}
