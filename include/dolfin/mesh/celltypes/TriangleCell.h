// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TRIANGLE_CELL_H
#define __DOLFIN_TRIANGLE_CELL_H

#include <dolfin/mesh/celltypes/CellType.h>
#include <dolfin/mesh/entities/Cell.h>

namespace dolfin
{

/**
 *  @class  TriangleCell
 *
 *  @brief  This class implements functionality for triangular meshes.
 *
 */

class TriangleCell : public CellType
{
  // UFC: Topological Dimension
  static size_t const TD = 2;

  // UFC: Number of Entities
  static size_t const NE[3][3];

  // UFC: Vertex Coordinates
  static real const VC[3][2];

  // UFC: Edge - Incident Vertices
  static size_t const EIV[3][2];

  // UFC: Edge - Non-Incident Vertices
  static size_t const ENV[3][1];

public:
  /// Specify cell type and facet type
  TriangleCell();

  ///
  ~TriangleCell() override = default;

  /// Clone pattern
  auto clone() const -> CellType * override
  {
    return new TriangleCell( *this );
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
  void create_entities( size_t ** e, size_t dim, size_t const * v ) const override;

  /// Order entities locally (connectivity 1-0, 2-0, 2-1)
  void order_entities( MeshTopology & topology, size_t i ) const override;

  /// Order vertices such that the facet is right-oriented w.r.t. facet normal
  void order_facet( size_t vertices[], Facet & facet ) const override;

  /// Return if mesh connectivities require ordering
  auto connectivity_needs_ordering( size_t d0, size_t d1 ) const -> bool override;

  /// Initialize mesh connectivities required by ordering
  void initialize_connectivities( Mesh & mesh ) const override;

  //--- REFINEMENT PATTERN ----------------------------------------------------

  /// Refine cell uniformly
  void
    refine_cell( Cell & cell, MeshEditor & editor, size_t & current_cell ) const override;

  /// Number of cells created by refinement pattern
  auto num_refined_cells() const -> size_t override;

  /// Number of vertices created by refinement pattern restricted to each
  /// entity of given topological dimensions
  auto num_refined_vertices( size_t dim ) const -> size_t override;

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
  void normal( Cell const & cell, size_t facet, real * n ) const override;

  /// Compute the area/length of given facet with respect to the cell
  auto facet_area( Cell const & cell, size_t facet ) const -> real override;

  /// Check if point p intersects the entity
  auto intersects( MeshEntity const & e, Point const & p ) const -> bool override;

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
inline auto TriangleCell::dim() const -> size_t
{
  return 2;
}

//-----------------------------------------------------------------------------
inline auto TriangleCell::num_entities( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= TD );
  return NE[2][dim];
}

//-----------------------------------------------------------------------------
inline auto TriangleCell::num_entities( size_t d0, size_t d1 ) const -> size_t
{
  dolfin_assert( d0 <= TD );
  dolfin_assert( d1 <= TD );
  return NE[d0][d1];
}

//-----------------------------------------------------------------------------
inline auto TriangleCell::num_vertices( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= TD );
  return NE[dim][0];
}

//-----------------------------------------------------------------------------
inline auto TriangleCell::orientation( Cell const & cell ) const -> size_t
{
  dolfin_assert( cell.type() == this->cell_type );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = cell.mesh().geometry();
  Array< size_t > const & vertices = cell.entities( 0 );
  real const *          v0       = geometry.x( vertices[0] );
  real const *          v1       = geometry.x( vertices[1] );
  real const *          v2       = geometry.x( vertices[2] );
  return ( ( ( v1[0] - v0[0] ) * ( v2[1] - v0[1] )
             - ( v1[1] - v0[1] ) * ( v2[0] - v0[0] ) )
               < 0.0
             ? 1
             : 0 );
}

//-----------------------------------------------------------------------------
inline auto TriangleCell::connectivity_needs_ordering( size_t d0, size_t d1 ) const -> bool
{
  dolfin_assert( d0 <= TD && d1 <= TD );
  return ( d0 > 0 && d0 > d1 );
}

//-----------------------------------------------------------------------------
inline void TriangleCell::initialize_connectivities( Mesh & mesh ) const
{
  mesh.init( 1, 0 );
  mesh.init( 2, 0 );
  mesh.init( 2, 1 );
}

//-----------------------------------------------------------------------------
inline auto TriangleCell::num_refined_cells() const -> size_t
{
  return 4;
}

//-----------------------------------------------------------------------------
inline auto TriangleCell::num_refined_vertices( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= TD );
  return ( dim > 1 ? 0 : 1 );
}

//-----------------------------------------------------------------------------
inline auto TriangleCell::volume( MeshEntity const & entity ) const -> real
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< size_t > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );

  switch ( geometry.dim() )
  {
    case 2:
      // Compute area of triangle embedded in R^2
      // Formula for volume from http://mathworld.wolfram.com
      return 0.5
             * std::abs( ( x0[0] * x1[1] + x0[1] * x2[0] + x1[0] * x2[1] )
                         - ( x2[0] * x1[1] + x2[1] * x0[0] + x1[0] * x0[1] ) );
      break;
    case 3:
      // Compute area of triangle embedded in R^3
      // Formula for volume from http://mathworld.wolfram.com
      return 0.5
             * std::sqrt(
               +std::pow( ( x0[1] * x1[2] + x0[2] * x2[1] + x1[1] * x2[2] )
                            - ( x2[1] * x1[2] + x2[2] * x0[1] + x1[1] * x0[2] ),
                          2 )
               + std::pow(
                 ( x0[2] * x1[0] + x0[0] * x2[2] + x1[2] * x2[0] )
                   - ( x2[2] * x1[0] + x2[0] * x0[2] + x1[2] * x0[0] ),
                 2 )
               + std::pow(
                 ( x0[0] * x1[1] + x0[1] * x2[0] + x1[0] * x2[1] )
                   - ( x2[0] * x1[1] + x2[1] * x0[0] + x1[0] * x0[1] ),
                 2 ) );
      break;
    default:
      error( "Volume of triangle only implemented for R^2 or R^3." );
      break;
  }
  return 0.0;
}

//-----------------------------------------------------------------------------
inline auto TriangleCell::diameter( MeshEntity const & entity ) const -> real
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< size_t > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );

  real e0 = 0.0;
  real e1 = 0.0;
  real e2 = 0.0;

  for ( size_t i = 0; i < geometry.dim(); ++i )
  {
    e0 += ( x1[i] - x0[i] ) * ( x1[i] - x0[i] );
    e1 += ( x2[i] - x1[i] ) * ( x2[i] - x1[i] );
    e2 += ( x0[i] - x2[i] ) * ( x0[i] - x2[i] );
  }

  return std::sqrt( std::max( std::max( e0, e1 ), e2 ) );
}

//-----------------------------------------------------------------------------
inline auto TriangleCell::circumradius( MeshEntity const & entity ) const -> real
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< size_t > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );

  real e0 = 0.0;
  real e1 = 0.0;
  real e2 = 0.0;

  for ( size_t i = 0; i < geometry.dim(); ++i )
  {
    e0 += ( x1[i] - x0[i] ) * ( x1[i] - x0[i] );
    e1 += ( x2[i] - x1[i] ) * ( x2[i] - x1[i] );
    e2 += ( x0[i] - x2[i] ) * ( x0[i] - x2[i] );
  }
  e0 = std::sqrt( e0 );
  e1 = std::sqrt( e1 );
  e2 = std::sqrt( e2 );

  //
  return 0.5
         * std::sqrt( ( e0 + e1 - e2 ) * ( e2 + e0 - e1 ) * ( e1 + e2 - e0 )
                      / ( e0 + e1 + e2 ) );
}

//-----------------------------------------------------------------------------
inline auto TriangleCell::inradius( MeshEntity const & entity ) const -> real
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< size_t > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );

  real e0 = 0.0;
  real e1 = 0.0;
  real e2 = 0.0;

  for ( size_t i = 0; i < geometry.dim(); ++i )
  {
    e0 += ( x1[i] - x0[i] ) * ( x1[i] - x0[i] );
    e1 += ( x2[i] - x1[i] ) * ( x2[i] - x1[i] );
    e2 += ( x0[i] - x2[i] ) * ( x0[i] - x2[i] );
  }
  e0 = std::sqrt( e0 );
  e1 = std::sqrt( e1 );
  e2 = std::sqrt( e2 );

  // Formula for circumradius from http://mathworld.wolfram.com
  // Using Heron's formula for the volume instead of calling volume()
  real const s = 0.5 * ( e0 + e1 + e2 );
  return 0.25 * e0 * e1 * e2
         / std::sqrt( s * ( s - e0 ) * ( s - e1 ) * ( s - e2 ) );
}

//-----------------------------------------------------------------------------
inline void TriangleCell::midpoint( MeshEntity const & entity, real * p ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< size_t > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );

  for ( size_t i = 0; i < geometry.dim(); ++i )
  {
    p[i] = ( x0[i] + x1[i] + x2[i] ) / 3.0;
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_TRIANGLE_CELL_H */
