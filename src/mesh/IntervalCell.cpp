// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-06-05
// Last changed: 2008-05-09
//
// Modified by Kristian Oelgaard, 2007.

#include <algorithm>
#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/IntervalCell.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
IntervalCell::IntervalCell() :
    CellType(interval, point)
{
}
//-----------------------------------------------------------------------------
uint IntervalCell::dim() const
{
  return 1;
}
//-----------------------------------------------------------------------------
uint IntervalCell::numEntities(uint dim) const
{
  switch (dim)
    {
    case 0:
      return 2;  // vertices
    case 1:
      return 1;  // cells
    default:
      error("Illegal topological dimension %d for interval.", dim);
      break;
    }

  return 0;
}
//-----------------------------------------------------------------------------
uint IntervalCell::numVertices(uint dim) const
{
  switch (dim)
    {
    case 0:
      return 1;  // vertices
    case 1:
      return 2;  // cells
    default:
      error("Illegal topological dimension %d for interval.", dim);
      break;
    }

  return 0;
}
//-----------------------------------------------------------------------------
uint IntervalCell::orientation(const Cell& cell) const
{
  Point v01 = Point(cell.entities(0)[1]) - Point(cell.entities(0)[0]);
  Point n(-v01.y(), v01.x());

  return (n.dot(v01) < 0.0 ? 1 : 0);
}
//-----------------------------------------------------------------------------
void IntervalCell::createEntities(uint** e, uint dim, uint const* v) const
{
  // We don't need to create any entities
  error("Don't know how to create entities of topological dimension %d.", dim);
}
//-----------------------------------------------------------------------------
void IntervalCell::orderEntities(Cell& cell) const
{
  // Sort i - j for i > j: 1 - 0

  // Get mesh topology
  MeshTopology& topology = cell.mesh().topology();

  // Sort local vertices in ascending order, connectivity 1 - 0
  if (topology(1, 0).size() > 0)
  {
    uint* cell_vertices = cell.entities(0);
    std::sort(cell_vertices, cell_vertices + 2);
  }
}
//-----------------------------------------------------------------------------
void IntervalCell::refineCell(Cell& cell, MeshEditor& editor,
                              uint& current_cell) const
{
  // Get vertices and edges
  uint const* v = cell.entities(0);
  uint const* e = cell.entities(1);
  dolfin_assert(v);
  dolfin_assert(e);

  // Get offset for new vertex indices
  uint const offset = cell.mesh().numVertices();

  // Compute indices for the three new vertices
  uint const e0 = offset + e[0];
  uint c0[2] = { v[0], e0 };
  uint c1[2] = { e0, v[1] };

  // Add the two new cells
  editor.addCell(current_cell++, &c0[0]);
  editor.addCell(current_cell++, &c1[0]);
}
//-----------------------------------------------------------------------------
real IntervalCell::volume(const MeshEntity& interval) const
{
  // Check that we get an interval
  dolfin_assert(interval.dim() == 1);

  // Get mesh geometry
  const MeshGeometry& geometry = interval.mesh().geometry();

  // Get the coordinates of the two vertices
  uint const* vertices = interval.entities(0);
  const real* x0 = geometry.x(vertices[0]);
  const real* x1 = geometry.x(vertices[1]);

  // Compute length of interval (line segment)
  real sum = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    const real dx = x1[i] - x0[i];
    sum += dx * dx;
  }

  return std::sqrt(sum);
}
//-----------------------------------------------------------------------------
real IntervalCell::diameter(const MeshEntity& interval) const
{
  // Diameter is same as volume for interval (line segment)
  return volume(interval);
}
//-----------------------------------------------------------------------------
real IntervalCell::circumradius(const MeshEntity& interval) const
{
  // Circumradius is same as volume for interval (line segment)
  return volume(interval);
}
//-----------------------------------------------------------------------------
real IntervalCell::normal(const Cell& cell, uint facet, uint i) const
{
  return normal(cell, facet)[i];
}
//-----------------------------------------------------------------------------
Point IntervalCell::normal(const Cell& cell, uint facet) const
{
  // Get mesh geometry
  const MeshGeometry& geometry = cell.mesh().geometry();

  // The normal vector is currently only defined for an interval in R^1
  if (geometry.dim() != 1) error(
      "The normal vector is only defined when the interval is in R^1");

  // Get the two vertices as points
  uint const* vertices = cell.entities(0);
  Point p0 = geometry.point(vertices[0]);
  Point p1 = geometry.point(vertices[1]);

  // Compute normal
  Point n = p0 - p1;
  if (facet == 1) n *= -1.0;

  // Normalize
  n /= n.norm();

  return n;
}
//-----------------------------------------------------------------------------
dolfin::real IntervalCell::facetArea(const Cell& cell, uint facet) const
{
  return 0.0;
}
//-----------------------------------------------------------------------------
bool IntervalCell::intersects(const MeshEntity& interval, const Point& p) const
{
  //FIXME: Due to constness inconsistency in Mesh
  Mesh * m = const_cast<Mesh *>(&interval.mesh());
  // Create points
  Point v0 = Vertex(*m, interval.entities(0)[0]).point();
  Point v1 = Vertex(*m, interval.entities(0)[1]).point();

  // Create vectors
  Point v01 = v1 - v0;
  Point vp0 = v0 - p;
  Point vp1 = v1 - p;

  // Check if the length of the sum of the two line segments vp0 and vp1 is
  // equal to the total length of the facet
  if (std::abs(v01.norm() - vp0.norm() - vp1.norm()) < DOLFIN_EPS)
  {
    return true;
  }
  else
  {
#if DEBUG
    message(2, "Point does not instersect with IntervalCell: "
            "epsilon = %f",
            std::abs(v01.norm() - vp0.norm() - vp1.norm()));
#endif
    return false;
  }
}
//-----------------------------------------------------------------------------
bool IntervalCell::intersects(const MeshEntity& interval, const Point& p1,
                              const Point& p2) const
{
  // FIXME: Not implemented
  error("Interval::intersects() not implemented");

  return false;
}
//-----------------------------------------------------------------------------
std::string IntervalCell::description() const
{
  std::string s = "interval (simplex of topological dimension 1)";
  return s;
}
//-----------------------------------------------------------------------------

}

