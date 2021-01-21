// Copyright (C) 2007-2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_POINT_CELL_H
#define __DOLFIN_POINT_CELL_H

#include <dolfin/math/basic.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshEditor.h>

namespace dolfin
{

/**
 *  @class  PointCell
 *
 *  @brief  This class implements functionality for point meshes.
 *
 */

class PointCell : public CellType
{
  // UFC: Topological Dimension
  static uint const TD = 0;

  // UFC: Number of Entities
  static uint const NE[1];

  // UFC: Vertex Coordinates
  static real const VC[1][1];

public:
  /// Specify cell type and facet type
  PointCell();

  ///
  ~PointCell() override = default;

  /// Clone pattern
  auto clone() const -> CellType * override
  {
    return new PointCell( *this );
  }

  /// Return topological dimension of cell
  auto dim() const -> uint override;

  /// Return number of entitites of given topological dimension
  auto num_entities( uint dim ) const -> uint override;

  /// Return number of entities of given topological dimensions
  auto num_entities( uint d0, uint d1 ) const -> uint override;

  /// Return number of vertices for entity of given topological dimension
  auto num_vertices( uint dim ) const -> uint override;

  /// Return orientation of the cell
  auto orientation( Cell const & cell ) const -> uint override;

  /// Create entities e of given topological dimension from vertices v
  void create_entities( uint ** e, uint dim, uint const * v ) const override;

  /// Order entities locally (connectivity 1-0, 2-0, 2-1)
  void order_entities( MeshTopology & topology, uint i ) const override;

  /// Order vertices such that the facet is right-oriented w.r.t. facet normal
  void order_facet( uint vertices[], Facet & facet ) const override;

  /// Return if mesh connectivities require ordering
  auto connectivity_needs_ordering( uint d0, uint d1 ) const -> bool override;

  /// Initialize mesh connectivities required by ordering
  void initialize_connectivities( Mesh & mesh ) const override;

  //--- REFINEMENT PATTERN ----------------------------------------------------

  /// Refine cell uniformly
  void refine_cell( Cell & cell, MeshEditor & editor,
                    uint & current_cell ) const override;

  /// Number of cells created by refinement pattern
  auto num_refined_cells() const -> uint override;

  /// Number of vertices created by refinement pattern restricted to each
  /// entity of given topological dimensions
  auto num_refined_vertices( uint dim ) const -> uint override;

  //---------------------------------------------------------------------------

  /// Compute (generalized) volume (area) of triangle
  auto volume( MeshEntity const & entity ) const -> real override;

  /// Compute diameter of triangle
  auto diameter( MeshEntity const & entity ) const -> real override;

  /// Compute circumradius of triangle
  auto circumradius( MeshEntity const & entity ) const -> real override;

  /// Compute inradius of interval
  auto inradius( MeshEntity const & entity ) const -> real override;

  /// Compute coordinates of midpoint
  void midpoint( MeshEntity const & entity, real * p ) const override;

  /// Compute of given facet with respect to the cell
  void normal( Cell const & cell, uint facet, real * n ) const override;

  /// Compute the area/length of given facet with respect to the cell
  auto facet_area( Cell const & cell, uint facet ) const -> real override;

  /// Check if point p intersects the cell
  auto intersects( MeshEntity const & e, Point const & p ) const -> bool override;

  /// Check if points line connecting p1 and p2 cuts the cell
  auto intersects( MeshEntity const & e,
                   Point const & p1, Point const & p2 ) const -> bool override;

  //--- REFERENCE CELL --------------------------------------------------------

  /// Create a mesh consisting of the reference cell
  void create_reference_cell( Mesh & mesh ) const override;

  /// Return coordinates of vertices in the reference cell
  auto reference_vertex( uint i ) const -> real const * override;

  //---------------------------------------------------------------------------

  /// Return description of cell type
  auto description() const -> std::string override;

  /// Display information
  void disp() const override;

  /// Check
  auto check( Cell & cell ) const -> bool override;

private:
};

//-----------------------------------------------------------------------------
inline auto PointCell::dim() const -> uint
{
  return 0;
}

//-----------------------------------------------------------------------------
inline auto PointCell::num_entities( uint dim ) const -> uint
{
  dolfin_assert( dim <= TD );
  MAYBE_UNUSED( dim );
  return 1;
}

//-----------------------------------------------------------------------------
inline auto PointCell::num_entities( uint d0, uint d1 ) const -> uint
{
  dolfin_assert( d0 <= TD );
  dolfin_assert( d1 <= TD );
  MAYBE_UNUSED( d0 );
  MAYBE_UNUSED( d1 );
  return 1;
}

//-----------------------------------------------------------------------------
inline auto PointCell::num_vertices( uint dim ) const -> uint
{
  dolfin_assert( dim <= TD );
  MAYBE_UNUSED( dim );
  return 1;
}

//-----------------------------------------------------------------------------
inline auto PointCell::orientation( Cell const & ) const -> uint
{
  return 0;
}

//-----------------------------------------------------------------------------
inline void
  PointCell::create_entities( uint ** e, uint dim, uint const * v ) const
{
  if ( dim > 0 )
  {
    error( "Invalid topological dimension for creation of entities: %d.", dim );
  }
  e[0][0] = v[0];
}

//-----------------------------------------------------------------------------
inline void PointCell::order_entities( MeshTopology &, uint ) const
{
  // do nothing
}

//-----------------------------------------------------------------------------
inline void PointCell::order_facet( uint[], Facet & ) const
{
  // Do nothing
}

//-----------------------------------------------------------------------------
inline auto PointCell::connectivity_needs_ordering( uint d0, uint d1 ) const -> bool
{
  dolfin_assert( d0 <= TD && d1 <= TD );
  MAYBE_UNUSED( d0 );
  MAYBE_UNUSED( d1 );
  return false;
}

//-----------------------------------------------------------------------------
inline void PointCell::initialize_connectivities( Mesh & ) const
{
  // Do nothing
}

//-----------------------------------------------------------------------------
inline void PointCell::refine_cell( Cell &       cell,
                                    MeshEditor & editor,
                                    uint &       current_cell ) const
{
  editor.add_cell( current_cell++, cell.entities( 0 ).data() );
}

//-----------------------------------------------------------------------------
inline auto PointCell::num_refined_cells() const -> uint
{
  return 1;
}

//-----------------------------------------------------------------------------
inline auto PointCell::num_refined_vertices( uint ) const -> uint
{
  return 1;
}

//-----------------------------------------------------------------------------
inline auto PointCell::volume( MeshEntity const & ) const -> real
{
  return 0.0;
}

//-----------------------------------------------------------------------------
inline auto PointCell::diameter( MeshEntity const & ) const -> real
{
  return 0.0;
}

//-----------------------------------------------------------------------------
inline auto PointCell::circumradius( MeshEntity const & ) const -> real
{
  return 0.0;
}

//-----------------------------------------------------------------------------
inline auto PointCell::inradius( MeshEntity const & ) const -> real
{
  return 0.0;
}

//-----------------------------------------------------------------------------
inline void PointCell::midpoint( MeshEntity const & entity, real * p ) const
{
  // Check that we get a point
  dolfin_assert( entity.dim() == 0 );
  dolfin_assert( entity.num_entities( 0 ) == 1 );
  real const * p0 = entity.mesh().geometry().x( entity.index() );
  std::copy( p0, p0 + entity.mesh().geometry_dimension(), p );
}

//-----------------------------------------------------------------------------
inline void PointCell::normal( Cell const &, uint, real * ) const
{
  error( "PointCell::normal() is undefined" );
}

//-----------------------------------------------------------------------------
inline auto PointCell::facet_area( Cell const &, uint ) const -> real
{
  return 0.0;
}

//-----------------------------------------------------------------------------
inline auto PointCell::intersects( MeshEntity const & e, Point const & p ) const -> bool
{
  return abscmp( p.dist( e.mesh().geometry().point( e.index() ) ), 0.0 );
}

//-----------------------------------------------------------------------------
inline auto PointCell::intersects( MeshEntity const &,
                                   Point const &,
                                   Point const & ) const -> bool
{
  error( "PointCell::intersects() not implemented." );
  return true;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_POINT_CELL_H */
