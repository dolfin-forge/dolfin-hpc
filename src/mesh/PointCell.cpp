// Copyright (C) 2007-2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Aurelien Larcher, 2015.
//
// First added:  2007-12-12
// Last changed: 2008-06-20

#include <dolfin/mesh/PointCell.h>

#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/Point.h>

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
    CellType(CellType::point, CellType::point)
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
  error("PointCell::orientation() undefined.");
  return 0;
}
//-----------------------------------------------------------------------------
void PointCell::create_entities(uint** e, uint dim, uint const* v) const
{
  error("PointCell::createEntities() undefined.");
}
//-----------------------------------------------------------------------------
void PointCell::order_entities(Cell& cell) const
{
  error("PointCell::orderEntities() undefined.");
}
//-----------------------------------------------------------------------------
void PointCell::refine_cell(Cell& cell, MeshEditor& editor,
                           uint& current_cell) const
{
  error("PointCell::refine_cell() undefined.");
}
//-----------------------------------------------------------------------------
uint PointCell::num_refined_cells() const
{
  error("PointCell::num_refined_cells() undefined.");
  return 0;
}
//-----------------------------------------------------------------------------
uint PointCell::num_refined_vertices(uint dim) const
{
  error("PointCell::num_refined_vertices() undefined.");
  return 0;
}
//-----------------------------------------------------------------------------
bool PointCell::refinement_needs_entities(uint dim) const
{
  error("PointCell::needs_entity_refined() undefined.");
  return false;
}
//-----------------------------------------------------------------------------
real PointCell::volume(MeshEntity const& entity) const
{
  error("PointCell::volume() undefined.");
  return 0.0;
}
//-----------------------------------------------------------------------------
real PointCell::diameter(MeshEntity const& entity) const
{
  error("PointCell::diameter() undefined.");
  return 0.0;
}
//-----------------------------------------------------------------------------
real PointCell::circumradius(MeshEntity const& entity) const
{
  error("PointCell::circumradius() undefined.");
  return 0.0;
}
//-----------------------------------------------------------------------------
Point PointCell::midpoint(MeshEntity const& entity) const
{
  // Check that we get a point
  dolfin_assert(entity.dim() == 0);
  dolfin_assert(entity.num_entities(0) == 1);

  return Point();
}
//-----------------------------------------------------------------------------
Point PointCell::normal(Cell const& cell, uint facet) const
{
  error("PointCell::normal() undefined.");
  return Point();
}
//-----------------------------------------------------------------------------
real PointCell::facet_area(Cell const& cell, uint facet) const
{
  error("PointCell::facetAread() undefined.");
  return 0.0;
}
//-----------------------------------------------------------------------------
bool PointCell::intersects(MeshEntity const& e, Point const& p) const
{
  error("PointCell::intersects() not implemented.");
  return true;
}
//-----------------------------------------------------------------------------
bool PointCell::intersects(MeshEntity const& e, Point const& p1,
                           Point const& p2) const
{
  error("PointCell::intersects() not implemented.");
  return true;
}
//-----------------------------------------------------------------------------
Mesh PointCell::create_reference_cell() const
{
  Mesh refcell;
  MeshEditor me(refcell, CellType::point, 1);
  me.init_vertices(1);
  me.add_vertex(0, VC[0]);
  me.init_cells(1);
  uint const cv0[1] = { 0 };
  me.add_cell(0, cv0);
  me.close();
  return refcell;
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
void PointCell::check(Cell& cell) const
{
  CellType::check(cell);
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
