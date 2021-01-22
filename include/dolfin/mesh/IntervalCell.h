// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_INTERVAL_CELL_H
#define __DOLFIN_INTERVAL_CELL_H

#include <dolfin/common/constants.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Vertex.h>

namespace dolfin
{

/**
 *  @class  IntervalCell
 *
 *  @brief  This class implements functionality for interval meshes.
 *
 */

class IntervalCell : public CellType
{
  // UFC: Topological Dimension
  static size_t const TD = 1;

  // UFC: Number of Entities
  static size_t const NE[2][2];

  // UFC: Vertex Coordinates
  static real const VC[2][1];

  // UFC: Edge - Incident Vertices
  static size_t const EIV[1][2];

public:
  /// Specify cell type and facet type
  IntervalCell();

  ///
  ~IntervalCell() override = default;

  /// Clone pattern
  auto clone() const -> CellType * override
  {
    return new IntervalCell( *this );
  }

  /// Return topological dimension of cell
  auto dim() const -> size_t override;

  /// Return number of entitites of given topological dimension
  auto num_entities( size_t dim ) const -> size_t override;

  /// Return number of entities of given topological dimensions
  auto num_entities( size_t d0, size_t d1 ) const -> size_t override;

  /// Return number of vertices for entity of given topological dimension
  auto num_vertices( size_t dim ) const -> size_t override;

  /// Return orientation of the cell
  auto orientation( Cell const & cell ) const -> size_t override;

  /// Create entities e of given topological dimension from vertices v
  void
    create_entities( size_t ** e, size_t dim, size_t const * v ) const override;

  /// Order entities locally (connectivity 1-0)
  void order_entities( MeshTopology & topology, size_t i ) const override;

  /// Order vertices such that the facet is right-oriented w.r.t. facet normal
  void order_facet( size_t vertices[], Facet & facet ) const override;

  /// Return if mesh connectivities require ordering
  auto connectivity_needs_ordering( size_t d0, size_t d1 ) const
    -> bool override;

  /// Initialize mesh connectivities required by ordering
  void initialize_connectivities( Mesh & mesh ) const override;

  //--- REFINEMENT PATTERN ----------------------------------------------------

  /// Refine cell uniformly
  void refine_cell( Cell &       cell,
                    MeshEditor & editor,
                    size_t &     current_cell ) const override;

  /// Number of cells created by refinement pattern
  auto num_refined_cells() const -> size_t override;

  /// Number of vertices created by refinement pattern restricted to each
  /// entity of given topological dimensions
  auto num_refined_vertices( size_t dim ) const -> size_t override;

  //---------------------------------------------------------------------------

  /// Compute (generalized) volume (length) of interval
  auto volume( MeshEntity const & entity ) const -> real override;

  /// Compute diameter of interval
  auto diameter( MeshEntity const & entity ) const -> real override;

  /// Compute circumradius of interval
  auto circumradius( MeshEntity const & entity ) const -> real override;

  /// Compute inradius of interval
  auto inradius( MeshEntity const & entity ) const -> real override;

  /// Compute coordinates of midpoint
  void midpoint( MeshEntity const & entity, real * p ) const override;

  /// Compute of given facet with respect to the cell
  void normal( Cell const & cell, size_t facet, real * n ) const override;

  /// Compute the area/length of given facet with respect to the cell
  auto facet_area( Cell const & cell, size_t facet ) const -> real override;

  /// Check if point p intersects the entity
  auto intersects( MeshEntity const & e, Point const & p ) const
    -> bool override;

  /// Check if points line connecting p1 and p2 cuts the entity
  auto intersects( MeshEntity const & e,
                   Point const &      p1,
                   Point const &      p2 ) const -> bool override;

  //--- REFERENCE CELL --------------------------------------------------------

  /// Create a mesh consisting of the reference cell
  void create_reference_cell( Mesh & mesh ) const override;

  /// Return coordinates of vertices in the reference cell
  auto reference_vertex( size_t i ) const -> real const * override;

  //---------------------------------------------------------------------------

  /// Return description of cell type
  auto description() const -> std::string override;

  /// Display information
  void disp() const override;

  /// Check
  auto check( Cell & cell ) const -> bool override;
};

//-----------------------------------------------------------------------------
inline auto IntervalCell::dim() const -> size_t
{
  return 1;
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::num_entities( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= TD );
  return NE[1][dim];
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::num_entities( size_t d0, size_t d1 ) const -> size_t
{
  dolfin_assert( d0 <= TD );
  dolfin_assert( d1 <= TD );
  return NE[d0][d1];
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::num_vertices( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= TD );
  return NE[dim][0];
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::orientation( Cell const & cell ) const -> size_t
{
  dolfin_assert( cell.type() == this->cell_type );
  Point v01 = Point( cell.entities( 0 )[1], 0.0, 0.0 )
              - Point( cell.entities( 0 )[0], 0.0, 0.0 );
  Point n( -v01[1], v01[0], 0.0 );

  return ( n.dot( v01 ) < 0.0 ? 1 : 0 );
}

//-----------------------------------------------------------------------------
inline void IntervalCell::order_facet( size_t[], Facet & ) const
{
  // Do nothing
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::connectivity_needs_ordering( size_t d0,
                                                       size_t d1 ) const -> bool
{
  dolfin_assert( d0 <= TD && d1 <= TD );
  return ( d0 == TD && d1 == 0 );
}

//-----------------------------------------------------------------------------
inline void IntervalCell::initialize_connectivities( Mesh & mesh ) const
{
  mesh.init( 1, 0 );
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::num_refined_cells() const -> size_t
{
  return 2;
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::num_refined_vertices( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= TD );
  MAYBE_UNUSED( dim );
  return 1;
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::volume( MeshEntity const & entity ) const -> real
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[1][0] );

  // Get mesh geometry
  MeshGeometry const & geometry = entity.mesh().geometry();

  // Get the coordinates of the two vertices
  Array< size_t > const & vertices = entity.entities( 0 );
  real const *            x0       = geometry.x( vertices[0] );
  real const *            x1       = geometry.x( vertices[1] );

  // Compute length of interval (line segment)
  real sum = 0.0;
  for ( size_t i = 0; i < geometry.dim(); ++i )
  {
    sum += ( x1[i] - x0[i] ) * ( x1[i] - x0[i] );
  }

  return std::sqrt( sum );
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::diameter( MeshEntity const & entity ) const -> real
{
  // Diameter is same as volume for interval (line segment)
  return volume( entity );
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::circumradius( MeshEntity const & entity ) const
  -> real
{
  // Circumradius is same as volume for interval (line segment)
  return volume( entity );
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::inradius( MeshEntity const & entity ) const -> real
{
  // Inradius is same as volume for interval (line segment)
  return volume( entity );
}

//-----------------------------------------------------------------------------
inline void IntervalCell::midpoint( MeshEntity const & entity, real * p ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[1][0] );

  MeshGeometry const &    geometry = entity.mesh().geometry();
  Array< size_t > const & vertices = entity.entities( 0 );
  real const *            x0       = geometry.x( vertices[0] );
  real const *            x1       = geometry.x( vertices[1] );
  size_t const            gdim     = geometry.dim();
  for ( size_t d = 0; d < gdim; ++d )
  {
    p[d] = 0.5 * ( x0[d] + x1[d] );
  }
}

//-----------------------------------------------------------------------------
inline void
  IntervalCell::normal( Cell const & cell, size_t facet, real * n ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  MeshGeometry const &    geometry = cell.mesh().geometry();
  Array< size_t > const & vertices = cell.entities( 0 );

  Point p0 = geometry.point( vertices[facet] );
  Point p1 = geometry.point( vertices[( facet + 1 ) % 2] );
  real  nn = p0.dist( p1 );

  for ( size_t d = 0; d < geometry.dim(); ++d )
  {
    n[d] /= nn;
  }
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::intersects( MeshEntity const & e,
                                      Point const &      p ) const -> bool
{
  dolfin_assert( e.dim() == TD );
  dolfin_assert( e.num_entities( 0 ) == NE[1][0] );

  // Get the coordinates of the vertices
  MeshGeometry const &    geometry = e.mesh().geometry();
  Array< size_t > const & vertices = e.entities( 0 );

  // Create points
  Point v0 = geometry.point( vertices[0] );
  Point v1 = geometry.point( vertices[1] );

  // Create vectors
  Point v01 = v1 - v0;
  Point vp0 = v0 - p;
  Point vp1 = v1 - p;

  // Check if the length of the sum of the two line segments vp0 and vp1 is
  // equal to the total length of the facet
  return ( std::abs( v01.norm() - vp0.norm() - vp1.norm() ) < DOLFIN_EPS );
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::facet_area( Cell const & cell, size_t ) const -> real
{
  dolfin_assert( cell.type() == this->cell_type );
  MAYBE_UNUSED( cell );
  return 0.0;
}

//-----------------------------------------------------------------------------
inline auto IntervalCell::reference_vertex( size_t i ) const -> real const *
{
  return &VC[i][0];
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_INTERVAL_CELL_H */
