// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/celltypes/IntervalCell.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/maybe_unused.h>
#include <dolfin/mesh/MeshEditor.h>

#include <algorithm>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

// UFC: Number of Entities
size_t const IntervalCell::NE[2][2] = { { 1, 0 }, { 2, 1 } };

// UFC: Vertex Coordinates
real const IntervalCell::VC[2][1] = { { 0.0 }, { 1.0 } };

// UFC: Edge - Incident Vertices
size_t const IntervalCell::EIV[1][2] = { { 0, 1 } };

//-----------------------------------------------------------------------------
IntervalCell::IntervalCell()
  : CellType( "interval", CellType::interval, CellType::point )
{
}
//-----------------------------------------------------------------------------
void IntervalCell::create_entities( size_t **      e,
                                    size_t         dim,
                                    size_t const * v ) const
{
  // We do not need to create any entities
  switch ( dim )
  {
    case 0: /* Create facets */
      e[0][0] = v[0];
      e[1][0] = v[1];
      break;
    default:
      error( "Invalid topological dimension for creation of entities: %d.",
             dim );
      break;
  }
}
//-----------------------------------------------------------------------------
void IntervalCell::order_entities( MeshTopology & topology, size_t i ) const
{
  // Sort i - j for i > j: 1 - 0
  dolfin_assert( topology.type().cellType() == this->cell_type );

  // Sort local vertices in ascending order, connectivity 1 - 0
  if ( topology.connectivity( 1, 0 ) )
  {
    Array< size_t > & cell_vertices = topology( 1, 0 )[i];
    std::sort( cell_vertices.data(), cell_vertices.data() + 2 );
  }
}
//-----------------------------------------------------------------------------
void IntervalCell::refine_cell( Cell &       cell,
                                MeshEditor & editor,
                                size_t &     current_cell ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  // Get vertices
  Array< size_t > const & v = cell.entities( 0 );
  dolfin_assert( !v.empty() );

  // Add midpoint vertex
  size_t const offset = cell.mesh().size( 0 );
  size_t const e0     = offset + cell.index();

  // Add the two new cells
  size_t const cv0[2] = { v[0], e0 };
  editor.add_cell( current_cell++, &cv0[0] );
  size_t const cv1[2] = { e0, v[1] };
  editor.add_cell( current_cell++, &cv1[0] );
}
//-----------------------------------------------------------------------------
auto IntervalCell::intersects( MeshEntity const & e,
                               Point const &,
                               Point const & ) const -> bool
{
  dolfin_assert( e.dim() == TD );
  dolfin_assert( e.num_entities( 0 ) == NE[1][0] );
  MAYBE_UNUSED( e );

  error( "Collision of interval with segment not implemented" );

  return false;
}
//-----------------------------------------------------------------------------
auto IntervalCell::description() const -> std::string
{
  return std::string( "interval (simplex of topological dimension 1)" );
}
//-----------------------------------------------------------------------------
void IntervalCell::create_reference_cell( Mesh & mesh ) const
{
  MeshEditor me( mesh, CellType::interval, 1, DOLFIN_COMM_SELF );
  me.init_vertices( 2 );
  me.add_vertex( 0, VC[0] );
  me.add_vertex( 1, VC[1] );
  me.init_cells( 1 );
  size_t const cv0[2] = { 0, 1 };
  me.add_cell( 0, cv0 );
  me.close();
}
//-----------------------------------------------------------------------------
void IntervalCell::disp() const
{
  message( "IntervalCell" );
  begin( "------------" );
  //---
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
auto IntervalCell::check( Cell & cell ) const -> bool
{
  bool ret = CellType::check( cell );

  // Check that cell vertices are in ascending order (so are edge vertices then)
  if ( cell.mesh().topology().connectivity( 1, 0 ) )
  {
    Array< size_t > const & cell_verts = cell.entities( 0 );
    dolfin_assert( not cell_verts.empty() );
    if ( cell_verts[1] < cell_verts[0] )
    {
      ret = false;
      warning( "Interval vertices are not in ascending order" );
    }
  }

  return ret;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
