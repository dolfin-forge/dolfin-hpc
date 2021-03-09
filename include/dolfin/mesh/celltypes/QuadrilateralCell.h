// Copyright (C) 2014 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_QUADRANGLE_CELL_H
#define __DOLFIN_QUADRANGLE_CELL_H

#include <dolfin/mesh/celltypes/CellType.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Facet.h>

namespace dolfin
{

/**
 *  @class  QuadrilateralCell
 *
 *  @brief  This class implements functionality for quadrilateral meshes.
 *
 */

class QuadrilateralCell : public CellType
{
  // UFC: Topological Dimension
  static size_t const TD = 2;

  // UFC: Number of Entities
  static size_t const NE[3][3];

  // UFC: Vertex Coordinates
  static real const VC[4][2];

  // UFC: Edge - Incident Vertices
  static size_t const EIV[4][2];

  // UFC: Edge - Non-Incident Vertices
  static size_t const ENV[4][2];

public:
  /// Specify cell type and facet type
  QuadrilateralCell();

  ///
  ~QuadrilateralCell() override = default;

  /// Clone pattern
  auto clone() const -> CellType * override
  {
    return new QuadrilateralCell( *this );
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

  /// Order entities locally (connectivity 1-0, 2-0, 2-1)
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

  /// Compute (generalized) volume (area) of quadrilateral
  auto volume( MeshEntity const & entity ) const -> real override;

  /// Compute diameter of quadrilateral
  auto diameter( MeshEntity const & entity ) const -> real override;

  /// Compute circumradius of quadrilateral
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

private:
  // Find local index of edge i according to ordering convention
  auto findEdge( size_t i, Cell const & cell ) const -> size_t;
};

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::dim() const -> size_t
{
  return 2;
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::num_entities( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= TD );
  return NE[2][dim];
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::num_entities( size_t d0, size_t d1 ) const
  -> size_t
{
  dolfin_assert( d0 <= TD );
  dolfin_assert( d1 <= TD );
  return NE[d0][d1];
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::num_vertices( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= TD );
  return NE[dim][0];
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::orientation( Cell const & cell ) const -> size_t
{
  dolfin_assert( cell.type() == this->cell_type );

  // Get the coordinates of vertices v0, v1 and v2
  MeshGeometry const &          geometry = cell.mesh().geometry();
  std::vector< size_t > const & vertices = cell.entities( 0 );
  real const *                  v0       = geometry.x( vertices[0] );
  real const *                  v1       = geometry.x( vertices[1] );
  real const *                  v2       = geometry.x( vertices[2] );

  // Check whether (v0v1, v0v2) is counter-clockwise
  return ( ( ( v1[0] - v0[0] ) * ( v2[1] - v0[1] )
             - ( v1[1] - v0[1] ) * ( v2[0] - v0[0] ) )
               < 0.0
             ? 1
             : 0 );
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::connectivity_needs_ordering( size_t d0,
                                                            size_t d1 ) const
  -> bool
{
  dolfin_assert( d0 <= TD && d1 <= TD );
  // Do not order cell - vertices connectivities
  return ( d0 > 0 && d0 > d1 ) && !( d0 == TD && d1 == 0 );
}

//-----------------------------------------------------------------------------
inline void QuadrilateralCell::initialize_connectivities( Mesh & mesh ) const
{
  mesh.init( 1, 0 );
  mesh.init( 2, 0 );
  mesh.init( 2, 1 );
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::num_refined_cells() const -> size_t
{
  return 4;
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::num_refined_vertices( size_t dim ) const
  -> size_t
{
  dolfin_assert( dim <= TD );
  MAYBE_UNUSED( dim );
  return 1;
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::volume( MeshEntity const & entity ) const -> real
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &          geometry = entity.mesh().geometry();
  std::vector< size_t > const & vertices = entity.entities( 0 );

  real const * a = geometry.x( vertices[0] );
  real const * b = geometry.x( vertices[1] );
  real const * c = geometry.x( vertices[2] );
  real const * d = geometry.x( vertices[3] );

  switch ( geometry.dim() )
  {
    case 2:
    {
      return 0.5
             * std::fabs( ( c[0] - a[0] ) * ( d[1] - b[1] )
                          - ( c[1] - a[1] ) * ( d[0] - b[0] ) );
    }
    break;
    case 3:
    {
      real c0 =
        ( c[1] - a[1] ) * ( d[2] - b[2] ) - ( c[2] - a[2] ) * ( d[1] - b[1] );
      real c1 =
        ( c[2] - a[2] ) * ( d[0] - b[0] ) - ( c[0] - a[0] ) * ( d[2] - b[2] );
      real c2 =
        ( c[0] - a[0] ) * ( d[1] - b[1] ) - ( c[1] - a[1] ) * ( d[0] - b[0] );
      return 0.5 * std::sqrt( c0 * c0 + c1 * c1 + c2 * c2 );
    }
    break;
  }

  return 0.0;
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::diameter( MeshEntity const & entity ) const
  -> real
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &          geometry = entity.mesh().geometry();
  std::vector< size_t > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );
  real const * x3 = geometry.x( vertices[3] );

  // Compute maximum diagonal length
  real d0 = 0.0;
  real d1 = 0.0;
  for ( size_t i = 0; i < geometry.dim(); ++i )
  {
    d0 += ( x3[i] - x0[i] ) * ( x3[i] - x0[i] );
    d1 += ( x2[i] - x1[i] ) * ( x2[i] - x1[i] );
  }
  return std::sqrt( std::max( d0, d1 ) );
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::circumradius( MeshEntity const & entity ) const
  -> real
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the four vertices
  MeshGeometry const &          geometry = entity.mesh().geometry();
  std::vector< size_t > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );
  real const * x3 = geometry.x( vertices[3] );

  // Compute circumradius assuming that the quadrilateral is cyclic
  real a = 0.0;
  real b = 0.0;
  real c = 0.0;
  real d = 0.0;
  for ( size_t i = 0; i < geometry.dim(); ++i )
  {
    a += ( x1[i] - x0[i] ) * ( x1[i] - x0[i] );
    b += ( x2[i] - x1[i] ) * ( x2[i] - x1[i] );
    c += ( x3[i] - x2[i] ) * ( x3[i] - x2[i] );
    d += ( x0[i] - x3[i] ) * ( x0[i] - x3[i] );
  }
  a = std::sqrt( a );
  b = std::sqrt( b );
  c = std::sqrt( c );
  d = std::sqrt( d );

  // Using Brahmagupta's formula for the volume
  real const s = 0.5 * ( a + b + c + d );
  return 0.25
         * std::sqrt( ( a * c + b * d ) * ( a * d + b * c ) * ( a * b + c * d )
                      / ( ( s - a ) * ( s - b ) * ( s - c ) * ( s - d ) ) );
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::inradius( MeshEntity const & entity ) const
  -> real
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the four vertices
  MeshGeometry const &          geometry = entity.mesh().geometry();
  std::vector< size_t > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );
  real const * x3 = geometry.x( vertices[3] );

  // Compute inradius assuming that the quadrilateral is tangent
  real a = 0.0;
  real b = 0.0;
  real c = 0.0;
  real d = 0.0;
  real p = 0.0;
  real q = 0.0;
  for ( size_t i = 0; i < geometry.dim(); ++i )
  {
    a += ( x1[i] - x0[i] ) * ( x1[i] - x0[i] );
    b += ( x2[i] - x1[i] ) * ( x2[i] - x1[i] );
    c += ( x3[i] - x2[i] ) * ( x3[i] - x2[i] );
    d += ( x0[i] - x3[i] ) * ( x0[i] - x3[i] );
    p += ( x2[i] - x0[i] ) * ( x2[i] - x0[i] );
    q += ( x3[i] - x1[i] ) * ( x3[i] - x1[i] );
  }
  a = std::sqrt( a );
  b = std::sqrt( b );
  c = std::sqrt( c );
  d = std::sqrt( d );
  p = std::sqrt( p );
  q = std::sqrt( q );

  return 0.5
         * std::sqrt( 4.0 * p * p * q * q
                      - ( a * a - b * b + c * c - d * d )
                          * ( a * a - b * b + c * c - d * d ) )
         / ( a + b + c + d );
}

//-----------------------------------------------------------------------------
inline void QuadrilateralCell::midpoint( MeshEntity const & entity,
                                         real *             p ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  MeshGeometry const &          geometry = entity.mesh().geometry();
  std::vector< size_t > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );
  real const * x3 = geometry.x( vertices[3] );

  for ( size_t d = 0; d < geometry.dim(); ++d )
  {
    p[d] = 0.25 * ( x0[d] + x1[d] + x2[d] + x3[d] );
  }
}

//-----------------------------------------------------------------------------
inline void
  QuadrilateralCell::normal( Cell const & cell, size_t facet, real * n ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  Cell &               c = const_cast< Cell & >( cell );
  Facet                f( c.mesh(), c.entities( 1 )[facet] );
  MeshGeometry const & geometry = cell.mesh().geometry();

  // Get coordinates of first non-incident vertex
  real const * p0 = geometry.x( c.entities( 0 )[ENV[facet][0]] );

  // Get coordinates of edge vertices
  std::vector< size_t > const & vertices = f.entities( 0 );

  real const * p1 = geometry.x( vertices[0] );
  real const * p2 = geometry.x( vertices[1] );

  switch ( geometry.dim() )
  {
    case 2:
    {
      n[0]          = p2[1] - p1[1];
      n[1]          = p1[0] - p2[0];
      real const nn = std::sqrt( n[0] * n[0] + n[1] * n[1] );
      n[0] /= nn;
      n[1] /= nn;
      // Flip direction of normal so it points outward
      if ( ( n[0] * ( p1[0] - p0[0] ) + n[1] * ( p1[1] - p0[1] ) ) < 0.0 )
      {
        n[0] *= -1.0;
        n[1] *= -1.0;
      }
    }
    break;
    case 3:
    {
      real ac = ( p2[0] - p1[0] ) * ( p2[0] - p0[0] )
                + ( p2[1] - p1[1] ) * ( p2[1] - p0[1] )
                + ( p2[2] - p1[2] ) * ( p2[2] - p0[2] );
      real ab = ( p2[0] - p1[0] ) * ( p1[0] - p0[0] )
                + ( p2[1] - p1[1] ) * ( p1[1] - p0[1] )
                + ( p2[2] - p1[2] ) * ( p1[2] - p0[2] );
      n[0]    = ( p1[0] - p0[0] ) * ac - ( p2[0] - p0[0] ) * ab;
      n[1]    = ( p1[1] - p0[1] ) * ac - ( p2[1] - p0[1] ) * ab;
      n[2]    = ( p1[2] - p0[2] ) * ac - ( p2[2] - p0[2] ) * ab;
      real nn = std::sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2] );
      n[0] /= nn;
      n[1] /= nn;
      n[2] /= nn;
      dolfin_assert( ( n[0] * ( p1[0] - p0[0] ) + n[1] * ( p1[1] - p0[1] )
                       + n[2] * ( p1[2] - p0[2] ) )
                     > 0.0 );
    }
    break;
    default:
      error( "QuadrilateralCell : facet normal not implemented in R^%u",
             geometry.dim() );
      break;
  }
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::facet_area( Cell const & cell,
                                           size_t       facet ) const -> real
{
  dolfin_assert( cell.type() == this->cell_type );

  Cell &                        c = const_cast< Cell & >( cell );
  Facet                         f( c.mesh(), c.entities( 1 )[facet] );
  MeshGeometry const &          geometry = cell.mesh().geometry();
  std::vector< size_t > const & vertices = f.entities( 0 );
  real const *                  p0       = geometry.x( vertices[0] );
  real const *                  p1       = geometry.x( vertices[1] );
  // Compute distance between vertices
  real         meas = 0.0;
  size_t const gdim = geometry.dim();
  for ( size_t d = 0; d < gdim; ++d )
  {
    meas += ( p1[d] - p0[d] ) * ( p1[d] - p0[d] );
  }
  return std::sqrt( meas );
}

//-----------------------------------------------------------------------------
inline auto QuadrilateralCell::reference_vertex( size_t i ) const
  -> real const *
{
  return &VC[i][0];
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_QUADRANGLE_CELL_H */
