// Copyright (C) 2007-2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/PointCell.h>

#include <dolfin/common/maybe_unused.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/Point.h>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
uint const PointCell::NE[1] =
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
