// Copyright (C) 2007-2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/celltypes/PointCell.h>

#include <dolfin/common/maybe_unused.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/entities/MeshEntity.h>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
size_t const PointCell::NE[1] = { 1 };

// UFC: Vertex Coordinates
real const PointCell::VC[1][1] = { { 0.0 } };

//-----------------------------------------------------------------------------
PointCell::PointCell()
  : CellType( "point", CellType::point, CellType::point )
{
}
//-----------------------------------------------------------------------------
void PointCell::create_reference_cell( Mesh & mesh ) const
{
  MeshEditor me( mesh, CellType::point, 1, DOLFIN_COMM_SELF );
  me.init_vertices( 1 );
  me.add_vertex( 0, VC[0] );
  me.init_cells( 1 );
  size_t const cv0[1] = { 0 };
  me.add_cell( 0, cv0 );
  me.close();
}
//-----------------------------------------------------------------------------
auto PointCell::reference_vertex( size_t i ) const -> real const *
{
  return &VC[i][0];
}
//-----------------------------------------------------------------------------
auto PointCell::description() const -> std::string
{
  return std::string( "point (simplex of topological dimension 0)" );
}
//-----------------------------------------------------------------------------
void PointCell::disp() const
{
  message( "PointCell" );
  begin( "---------" );
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
auto PointCell::check( Cell & cell ) const -> bool
{
  return CellType::check( cell );
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
