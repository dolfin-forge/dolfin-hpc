// Copyright (C) 2007-2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Aurelien Larcher, 2015.
//
// First added:  2007-12-12
// Last changed: 2008-06-20

#include <dolfin/mesh/PointCell.h>

#include <dolfin/math/basic.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshEditor.h>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
uint const PointCell::NE[1] =
{ 1 };

// UFC: Number of Vertices (per entity)
uint const PointCell::NV[1] =
{ 1 };

// UFC: Vertex Coordinates
real const PointCell::VC[1][1] =
{ { 0.0 } };

//-----------------------------------------------------------------------------
PointCell::PointCell() :
    CellType("point", CellType::point, CellType::point)
{
}
//-----------------------------------------------------------------------------
PointCell::~PointCell()
{
}
//-----------------------------------------------------------------------------
uint PointCell::dim() const
{
  return 0;
}
//-----------------------------------------------------------------------------
uint PointCell::num_entities(uint dim) const
{
  dolfin_assert(dim <= TD);
  return 1;
}
//-----------------------------------------------------------------------------
uint PointCell::num_vertices(uint dim) const
{
  dolfin_assert(dim <= TD);
  return 1;
}
//-----------------------------------------------------------------------------
uint PointCell::orientation(Cell const& cell) const
{
  return 0;
}
//-----------------------------------------------------------------------------
void PointCell::create_entities(uint** e, uint dim, uint const* v) const
{
  if(dim > 0)
  {
    error("Invalid topological dimension for creation of entities: %d.", dim);
  }
  e[0][0] = v[0];
}
//-----------------------------------------------------------------------------
void PointCell::order_entities(MeshTopology& topology, uint i) const
{
  // do nothing
}
//-----------------------------------------------------------------------------
void PointCell::order_facet(uint vertices[], Facet& facet) const
{
  // Do nothing
}
//-----------------------------------------------------------------------------
bool PointCell::connectivity_needs_ordering(uint d0, uint d1) const
{
  dolfin_assert(d0 <= TD && d1 <= TD);
  return false;
}
//-----------------------------------------------------------------------------
void PointCell::initialize_connectivities(Mesh& mesh) const
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void PointCell::refine_cell(Cell& cell, MeshEditor& editor,
                           uint& current_cell) const
{
  editor.add_cell(current_cell++, cell.entities(0));
}
//-----------------------------------------------------------------------------
uint PointCell::num_refined_cells() const
{
  return 1;
}
//-----------------------------------------------------------------------------
uint PointCell::num_refined_vertices(uint dim) const
{
  return 1;
}
//-----------------------------------------------------------------------------
real PointCell::volume(MeshEntity const& entity) const
{
  return 0.0;
}
//-----------------------------------------------------------------------------
real PointCell::diameter(MeshEntity const& entity) const
{
  return 0.0;
}
//-----------------------------------------------------------------------------
real PointCell::circumradius(MeshEntity const& entity) const
{
  return 0.0;
}
//-----------------------------------------------------------------------------
real PointCell::inradius(MeshEntity const& entity) const
{
  return 0.0;
}
//-----------------------------------------------------------------------------
void PointCell::midpoint(MeshEntity const& entity, real * p) const
{
  // Check that we get a point
  dolfin_assert(entity.dim() == 0);
  dolfin_assert(entity.num_entities(0) == 1);
  real const * p0 = entity.mesh().geometry().x(entity.index());
  std::copy(p0, p0 + entity.mesh().geometry_dimension(), p);
}
//-----------------------------------------------------------------------------
void PointCell::normal(Cell const& cell, uint facet, real * n) const
{
  error("PointCell::normal() is undefined");
}
//-----------------------------------------------------------------------------
real PointCell::facet_area(Cell const& cell, uint facet) const
{
  return 0.0;
}
//-----------------------------------------------------------------------------
bool PointCell::intersects(MeshEntity const& e, Point const& p) const
{
  return abscmp(p.dist(e.mesh().geometry().point(e.index())), 0.0);
}
//-----------------------------------------------------------------------------
bool PointCell::intersects(MeshEntity const& e, Point const& p1,
                           Point const& p2) const
{
  error("PointCell::intersects() not implemented.");
  return true;
}
//-----------------------------------------------------------------------------
void PointCell::create_reference_cell(Mesh& mesh) const
{
  MeshEditor me(mesh, CellType::point, 1, DOLFIN_COMM_SELF);
  me.init_vertices(1);
  me.add_vertex(0, VC[0]);
  me.init_cells(1);
  uint const cv0[1] = { 0 };
  me.add_cell(0, cv0);
  me.close();
}
//-----------------------------------------------------------------------------
real const * PointCell::reference_vertex(uint i) const
{
  return &VC[i][0];
}
//-----------------------------------------------------------------------------
std::string PointCell::description() const
{
  return std::string("point (simplex of topological dimension 0)");
}
//-----------------------------------------------------------------------------
void PointCell::disp() const
{
  message("PointCell");
  begin(  "---------");
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
bool PointCell::check(Cell& cell) const
{
  return CellType::check(cell);
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
