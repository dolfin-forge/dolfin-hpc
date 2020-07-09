// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TRIANGLE_CELL_H
#define __DOLFIN_TRIANGLE_CELL_H

#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Cell.h>

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
  static uint const TD = 2;

  // UFC: Number of Entities
  static uint const NE[3][3];

  // UFC: Vertex Coordinates
  static real const VC[3][2];

  // UFC: Edge - Incident Vertices
  static uint const EIV[3][2];

  // UFC: Edge - Non-Incident Vertices
  static uint const ENV[3][1];

public:
  /// Specify cell type and facet type
  TriangleCell();

  ///
  ~TriangleCell();

  /// Clone pattern
  CellType * clone() const
  {
    return new TriangleCell( *this );
  }

  /// Return topological dimension of cell
  uint dim() const;

  /// Return number of entitites of given topological dimension
  uint num_entities( uint dim ) const;

  /// Return number of entities of given topological dimensions
  uint num_entities( uint d0, uint d1 ) const;

  /// Return number of vertices for entity of given topological dimension
  uint num_vertices( uint dim ) const;

  /// Return orientation of the cell
  uint orientation( Cell const & cell ) const;

  /// Create entities e of given topological dimension from vertices v
  void create_entities( uint ** e, uint dim, uint const * v ) const;

  /// Order entities locally (connectivity 1-0, 2-0, 2-1)
  void order_entities( MeshTopology & topology, uint i ) const;

  /// Order vertices such that the facet is right-oriented w.r.t. facet normal
  void order_facet( uint vertices[], Facet & facet ) const;

  /// Return if mesh connectivities require ordering
  bool connectivity_needs_ordering( uint d0, uint d1 ) const;

  /// Initialize mesh connectivities required by ordering
  void initialize_connectivities( Mesh & mesh ) const;

  //--- REFINEMENT PATTERN ----------------------------------------------------

  /// Refine cell uniformly
  void
    refine_cell( Cell & cell, MeshEditor & editor, uint & current_cell ) const;

  /// Number of cells created by refinement pattern
  uint num_refined_cells() const;

  /// Number of vertices created by refinement pattern restricted to each
  /// entity of given topological dimensions
  uint num_refined_vertices( uint dim ) const;

  //---------------------------------------------------------------------------

  /// Compute (generalized) volume (area) of triangle
  real volume( MeshEntity const & entity ) const;

  /// Compute diameter of triangle
  real diameter( MeshEntity const & entity ) const;

  /// Compute circumradius of triangle
  real circumradius( MeshEntity const & entity ) const;

  /// Compute inradius of interval
  real inradius( MeshEntity const & entity ) const;

  /// Compute coordinates of midpoint
  void midpoint( MeshEntity const & entity, real * p ) const;

  /// Compute of given facet with respect to the cell
  void normal( Cell const & cell, uint facet, real * n ) const;

  /// Compute the area/length of given facet with respect to the cell
  real facet_area( Cell const & cell, uint facet ) const;

  /// Check if point p intersects the entity
  bool intersects( MeshEntity const & e, Point const & p ) const;

  /// Check if points line connecting p1 and p2 cuts the entity
  bool intersects( MeshEntity const & e,
                   Point const &      p1,
                   Point const &      p2 ) const;

  //--- REFERENCE CELL --------------------------------------------------------

  /// Create a mesh consisting of the reference cell
  void create_reference_cell( Mesh & mesh ) const;

  /// Return coordinates of vertices in the reference cell
  real const * reference_vertex( uint i ) const;

  //---------------------------------------------------------------------------

  /// Return description of cell type
  std::string description() const;

  /// Display information
  void disp() const;

  /// Check
  bool check( Cell & cell ) const;

private:
  // Find local index of edge i according to ordering convention
  uint findEdge( uint i, Cell const & cell ) const;
};

//-----------------------------------------------------------------------------
inline uint TriangleCell::dim() const
{
  return 2;
}

//-----------------------------------------------------------------------------
inline uint TriangleCell::num_entities( uint dim ) const
{
  dolfin_assert( dim <= TD );
  return NE[2][dim];
}

//-----------------------------------------------------------------------------
inline uint TriangleCell::num_entities( uint d0, uint d1 ) const
{
  dolfin_assert( d0 <= TD );
  dolfin_assert( d1 <= TD );
  return NE[d0][d1];
}

//-----------------------------------------------------------------------------
inline uint TriangleCell::num_vertices( uint dim ) const
{
  dolfin_assert( dim <= TD );
  return NE[dim][0];
}

//-----------------------------------------------------------------------------
inline uint TriangleCell::orientation( Cell const & cell ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = cell.mesh().geometry();
  Array< uint > const & vertices = cell.entities( 0 );
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
inline bool TriangleCell::connectivity_needs_ordering( uint d0, uint d1 ) const
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
inline uint TriangleCell::num_refined_cells() const
{
  return 4;
}

//-----------------------------------------------------------------------------
inline uint TriangleCell::num_refined_vertices( uint dim ) const
{
  dolfin_assert( dim <= TD );
  return ( dim > 1 ? 0 : 1 );
}

//-----------------------------------------------------------------------------
inline real TriangleCell::volume( MeshEntity const & entity ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities( 0 );

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
inline real TriangleCell::diameter( MeshEntity const & entity ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );

  real e0 = 0.0;
  real e1 = 0.0;
  real e2 = 0.0;

  for ( uint i = 0; i < geometry.dim(); ++i )
  {
    e0 += ( x1[i] - x0[i] ) * ( x1[i] - x0[i] );
    e1 += ( x2[i] - x1[i] ) * ( x2[i] - x1[i] );
    e2 += ( x0[i] - x2[i] ) * ( x0[i] - x2[i] );
  }

  return std::sqrt( std::max( std::max( e0, e1 ), e2 ) );
}

//-----------------------------------------------------------------------------
inline real TriangleCell::circumradius( MeshEntity const & entity ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );

  real e0 = 0.0;
  real e1 = 0.0;
  real e2 = 0.0;

  for ( uint i = 0; i < geometry.dim(); ++i )
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
inline real TriangleCell::inradius( MeshEntity const & entity ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[2][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );

  real e0 = 0.0;
  real e1 = 0.0;
  real e2 = 0.0;

  for ( uint i = 0; i < geometry.dim(); ++i )
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
  Array< uint > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );

  for ( uint i = 0; i < geometry.dim(); ++i )
  {
    p[i] = ( x0[i] + x1[i] + x2[i] ) / 3.0;
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_TRIANGLE_CELL_H */
