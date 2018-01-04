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
#include <dolfin/mesh/MeshEditor.h>

#include <algorithm>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
uint const IntervalCell::NE[2][2] =
{ { 1, 0 }, { 2, 1 } };

// UFC: Vertex Coordinates
real const IntervalCell::VC[2][1] =
{ { 0.0 }, { 1.0 } };

// UFC: Edge - Incident Vertices
uint const IntervalCell::EIV[1][2] =
{ { 0, 1 } };

//-----------------------------------------------------------------------------
IntervalCell::IntervalCell() :
    CellType("interval", CellType::interval, CellType::point)
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
  return NE[1][dim];
}
//-----------------------------------------------------------------------------
uint IntervalCell::num_entities(uint d0, uint d1) const
{
  dolfin_assert(d0 <= TD);
  dolfin_assert(d1 <= TD);
  return NE[d0][d1];
}
//-----------------------------------------------------------------------------
uint IntervalCell::num_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return NE[dim][0];
}
//-----------------------------------------------------------------------------
uint IntervalCell::orientation(Cell const& cell) const
{
  dolfin_assert(cell.type() == this->cell_type);
  Point v01 = Point(cell.entities(0)[1]) - Point(cell.entities(0)[0]);
  Point n(-v01[1], v01[0]);

  return (n.dot(v01) < 0.0 ? 1 : 0);
}
//-----------------------------------------------------------------------------
void IntervalCell::create_entities(uint** e, uint dim, uint const* v) const
{
  // We do not need to create any entities
  switch(dim)
  {
    case 0: /* Create facets */
      e[0][0] = v[0];
      e[1][0] = v[1];
      break;
    default:
      error("Invalid topological dimension for creation of entities: %d.", dim);
      break;
  }
}
//-----------------------------------------------------------------------------
void IntervalCell::order_entities(MeshTopology& topology, uint i) const
{
  // Sort i - j for i > j: 1 - 0
  dolfin_assert(topology.type(i).cellType() == this->cell_type);

  // Sort local vertices in ascending order, connectivity 1 - 0
  if (topology.is_computed(1, 0))
  {
    uint* cell_vertices = topology(1, 0)(i);
    std::sort(cell_vertices, cell_vertices + 2);
  }
}
//-----------------------------------------------------------------------------
void IntervalCell::order_facet(uint vertices[], Facet& facet) const
{
  // Do nothing
}
//-----------------------------------------------------------------------------
bool IntervalCell::connectivity_needs_ordering(uint d0, uint d1) const
{
  dolfin_assert(d0 <= TD && d1 <= TD);
  return (d0 == TD && d1 == 0);
}
//-----------------------------------------------------------------------------
void IntervalCell::initialize_connectivities(Mesh& mesh) const
{
  mesh.init(1, 0);
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
  uint const offset = cell.mesh().size(0);
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
real IntervalCell::volume(MeshEntity const& entity) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[1][0]);

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
real IntervalCell::inradius(MeshEntity const& entity) const
{
  // Inradius is same as volume for interval (line segment)
  return volume(entity);
}
//-----------------------------------------------------------------------------
void IntervalCell::midpoint(MeshEntity const& entity, real * p) const
{
  dolfin_assert(entity.dim() == TD);
  dolfin_assert(entity.num_entities(0) == NE[1][0]);

  MeshGeometry const& geometry = entity.mesh().geometry();
  uint const * vertices = entity.entities(0);
  real const * x0 = geometry.x(vertices[0]);
  real const * x1 = geometry.x(vertices[1]);
  uint const gdim = geometry.dim();
  for (uint d = 0; d < gdim; ++d)
  {
    p[d] = 0.5 * ( x0[d] + x1[d] );
  }
}
//-----------------------------------------------------------------------------
void IntervalCell::normal(Cell const& cell, uint facet, real * n) const
{
  dolfin_assert(cell.type() == this->cell_type);

  MeshGeometry const& geometry = cell.mesh().geometry();
  uint const * vertices = cell.entities(0);
  real const * p0 = geometry.x(vertices[facet]);
  real const * p1 = geometry.x(vertices[(facet + 1) % 2]);
  uint const gdim = geometry.dim();
  real nn = 0.0;
  for (uint d = 0; d < gdim; ++d)
  {
    n[d] = p0[d] - p1[d];
    nn += n[d] * n[d];
  }
  nn = std::sqrt(nn);
  for (uint d = 0; d < gdim; ++d)
  {
    n[d] /= nn;
  }
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
  dolfin_assert(e.num_entities(0) == NE[1][0]);

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
  dolfin_assert(e.num_entities(0) == NE[1][0]);

  error("Collision of interval with segment not implemented");

  return false;
}
//-----------------------------------------------------------------------------
std::string IntervalCell::description() const
{
  return std::string("interval (simplex of topological dimension 1)");
}
//-----------------------------------------------------------------------------
void IntervalCell::create_reference_cell(Mesh& mesh) const
{
  MeshEditor me(mesh, CellType::interval, 1, DOLFIN_COMM_SELF);
  me.init_vertices(2);
  me.add_vertex(0, VC[0]);
  me.add_vertex(1, VC[1]);
  me.init_cells(1);
  uint const cv0[2] = { 0, 1 };
  me.add_cell(0, cv0);
  me.close();
}
//-----------------------------------------------------------------------------
real const * IntervalCell::reference_vertex(uint i) const
{
  return &VC[i][0];
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
bool IntervalCell::check(Cell& cell) const
{
  bool ret = CellType::check(cell);

  // Check that cell vertices are in ascending order (so are edge vertices then)
  if (cell.mesh().topology().is_computed(1, 0))
  {
    uint* cell_verts = cell.entities(0);
    dolfin_assert(cell_verts);
    if (cell_verts[1] < cell_verts[0])
    {
      ret = false;
      warning("Interval vertices are not in ascending order");
    }
  }

  return ret;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
