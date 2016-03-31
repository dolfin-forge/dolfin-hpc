// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-06-05
// Last changed: 2008-05-09
//
// Modified by Kristian Oelgaard, 2007.
// Modified by Aurelien Larcher, 2015.
//

#include <dolfin/mesh/IntervalCell.h>

#include <dolfin/common/constants.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

#include <algorithm>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
uint const IntervalCell::NE[2] =
{ 2, 1 };

// UFC: Number of Vertices (per entity)
uint const IntervalCell::NV[2] =
{ 1, 2 };

// UFC: Vertex Coordinates
real const IntervalCell::VC[2][1] =
{ { 0.0 }, { 1.0 } };

// UFC: Edge - Incident Vertices
uint const IntervalCell::EIV[1][2] =
{ { 0, 1 } };

//-----------------------------------------------------------------------------
IntervalCell::IntervalCell() :
    CellType(CellType::interval, CellType::point)
{
}
//-----------------------------------------------------------------------------
IntervalCell::~IntervalCell()
{
}
//-----------------------------------------------------------------------------
uint IntervalCell::dim() const
{
  return 1;
}
//-----------------------------------------------------------------------------
uint IntervalCell::num_entities(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NE[dim];
}
//-----------------------------------------------------------------------------
uint IntervalCell::num_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NV[dim];
}
//-----------------------------------------------------------------------------
uint IntervalCell::orientation(Cell const& cell) const
{
  dolfin_assert(cell.type() == this->cell_type);
  Point v01 = Point(cell.entities(0)[1]) - Point(cell.entities(0)[0]);
  Point n(-v01.y(), v01.x());

  return (n.dot(v01) < 0.0 ? 1 : 0);
}
//-----------------------------------------------------------------------------
void IntervalCell::create_entities(uint** e, uint dim, uint const* v) const
{
  // We do not need to create any entities
  error("Invalid topological dimension for creation of entities: %d.", dim);
}
//-----------------------------------------------------------------------------
void IntervalCell::order_entities(Cell& cell) const
{
  // Sort i - j for i > j: 1 - 0
  dolfin_assert(cell.type() == this->cell_type);

  // Get mesh topology
  MeshTopology const& topology = cell.mesh().topology();

  // Sort local vertices in ascending order, connectivity 1 - 0
  if (topology(1, 0).size() > 0)
  {
    uint* cell_vertices = cell.entities(0);
    std::sort(cell_vertices, cell_vertices + 2);
  }
}
//-----------------------------------------------------------------------------
void IntervalCell::refine_cell(Cell& cell, MeshEditor& editor,
                              uint& current_cell) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get vertices
  uint const* v = cell.entities(0);
  dolfin_assert(v);

  // Add midpoint vertex
  uint const offset = cell.mesh().numVertices();
  uint const e0 = offset + cell.index();

  // Add the two new cells
  uint const cv0[2] = { v[0], e0 };
  editor.add_cell(current_cell++, &cv0[0]);
  uint const cv1[2] = { e0, v[1] };
  editor.add_cell(current_cell++, &cv1[0]);
}
//-----------------------------------------------------------------------------
uint IntervalCell::num_refined_cells() const
{
  return 2;
}
//-----------------------------------------------------------------------------
uint IntervalCell::num_refined_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return 1;
}
//-----------------------------------------------------------------------------
bool IntervalCell::refinement_needs_entities(uint dim) const
{
  dolfin_assert(dim <= TD);
  return true;
}
//-----------------------------------------------------------------------------
real IntervalCell::volume(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.numEntities(0) == NE[0]);

  // Get mesh geometry
  MeshGeometry const& geometry = entity.mesh().geometry();

  // Get the coordinates of the two vertices
  uint const* vertices = entity.entities(0);
  real const* x0 = geometry.x(vertices[0]);
  real const* x1 = geometry.x(vertices[1]);

  // Compute length of interval (line segment)
  real sum = 0.0;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    sum += (x1[i] - x0[i]) * (x1[i] - x0[i]);
  }

  return std::sqrt(sum);
}
//-----------------------------------------------------------------------------
real IntervalCell::diameter(MeshEntity const& entity) const
{
  // Diameter is same as volume for interval (line segment)
  return volume(entity);
}
//-----------------------------------------------------------------------------
real IntervalCell::circumradius(MeshEntity const& entity) const
{
  // Circumradius is same as volume for interval (line segment)
  return volume(entity);
}
//-----------------------------------------------------------------------------
Point IntervalCell::midpoint(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.numEntities(0) == NE[0]);

  // Get the coordinates of the vertices
  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const* vertices = entity.entities(0);
  real const* x0 = geometry.x(vertices[0]);
  real const* x1 = geometry.x(vertices[1]);
  Point p;
  for (uint i = 0; i < geometry.dim(); ++i)
  {
    p[i] = 0.5 * ( x0[i] + x1[i] );
  }
  return p;
}
//-----------------------------------------------------------------------------
Point IntervalCell::normal(Cell const& cell, uint facet) const
{
  dolfin_assert(cell.type() == this->cell_type);

  // Get mesh geometry
  MeshGeometry const& geometry = cell.mesh().geometry();
  if (geometry.dim() != 1)
  {
    error("The normal vector is only defined when the interval is in R^1");
  }

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
real IntervalCell::facet_area(Cell const& cell, uint facet) const
{
  dolfin_assert(cell.type() == this->cell_type);
  return 0.0;
}
//-----------------------------------------------------------------------------
bool IntervalCell::intersects(MeshEntity const& e, Point const& p) const
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.numEntities(0) == NE[0]);

  // Get the coordinates of the vertices
  MeshGeometry const& geometry = e.mesh().geometry();
  uint const* vertices = e.entities(0);
  real const* x0 = geometry.x(vertices[0]);
  real const* x1 = geometry.x(vertices[1]);

  // Create points
  Point v0;
  std::memcpy(&v0[0], x0, geometry.dim()*sizeof(real));
  Point v1;
  std::memcpy(&v1[0], x1, geometry.dim()*sizeof(real));

  // Create vectors
  Point v01 = v1 - v0;
  Point vp0 = v0 - p;
  Point vp1 = v1 - p;

  // Check if the length of the sum of the two line segments vp0 and vp1 is
  // equal to the total length of the facet
  return ( std::abs(v01.norm() - vp0.norm() - vp1.norm()) < DOLFIN_EPS );
}
//-----------------------------------------------------------------------------
bool IntervalCell::intersects(MeshEntity const& e, Point const& p1,
                              Point const& p2) const
{
  dolfin_assert(e.dim() == TD);
  dolfin_assert(e.numEntities(0) == NE[0]);

  error("Collision of interval with segment not implemented");

  return false;
}
//-----------------------------------------------------------------------------
std::string IntervalCell::description() const
{
  return std::string("interval (simplex of topological dimension 1)");
}
//-----------------------------------------------------------------------------
Mesh IntervalCell::create_reference_cell() const
{
  Mesh refcell;
  MeshEditor me(refcell, CellType::interval, 1);
  me.init_vertices(2);
  me.add_vertex(0, VC[0]);
  me.add_vertex(1, VC[1]);
  me.init_cells(1);
  uint const cv0[2] = { 0, 1 };
  me.add_cell(0, cv0);
  me.close();
  return refcell;
}
//-----------------------------------------------------------------------------
void IntervalCell::disp() const
{
  message("IntervalCell");
  begin(  "------------");
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
void IntervalCell::check(Cell& cell) const
{
  CellType::check(cell);
  // Check that cell vertices are in ascending order (so are edge vertices then)
  uint* cell_verts = cell.entities(0);
  dolfin_assert(cell_verts != NULL);
  if (cell_verts[1] < cell_verts[0])
  {
    error("Interval vertices are not in ascending order");
  }
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
