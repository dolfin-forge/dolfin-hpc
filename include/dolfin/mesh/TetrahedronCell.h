// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TETRAHEDRON_CELL_H
#define __DOLFIN_TETRAHEDRON_CELL_H

#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Cell.h>

namespace dolfin
{

/**
 *  @class  TetrahedronCell
 *
 *  @brief  This class implements functionality for tetrahedral meshes.
 *
 */

class TetrahedronCell : public CellType
{
  // UFC: Topological Dimension
  static uint const TD = 3;

  // UFC: Number of Entities
  static uint const NE[4][4];

  // UFC: Vertex Coordinates
  static real const VC[4][3];

  // UFC: Edge - Incident Vertices
  static uint const EIV[6][2];

  // UFC: Edge - Non-Incident Vertices
  static uint const ENV[6][2];

  // UFC: Face - Incident Vertices
  static uint const FIV[4][3];

  // UFC: Face - Non-Incident Vertices
  static uint const FNV[4][1];

public:
  /// Specify cell type and facet type
  TetrahedronCell();

  ///
  ~TetrahedronCell() override = default;

  /// Clone pattern
  CellType * clone() const override
  {
    return new TetrahedronCell( *this );
  }

  /// Return topological dimension of cell
  uint dim() const override;

  /// Return number of entitites of given topological dimension
  uint num_entities( uint dim ) const override;

  /// Return number of entities of given topological dimensions
  uint num_entities( uint d0, uint d1 ) const override;

  /// Return number of vertices for entity of given topological dimension
  uint num_vertices( uint dim ) const override;

  /// Return orientation of the cell
  uint orientation( Cell const & cell ) const override;

  /// Create entities e of given topological dimension from vertices v
  void create_entities( uint ** e, uint dim, uint const * v ) const override;

  /// Order entities locally (connectivity 1-0, 2-0, 2-1, 3-0, 3-1, 3-2)
  void order_entities( MeshTopology & topology, uint i ) const override;

  /// Order vertices such that the facet is right-oriented w.r.t. facet normal
  void order_facet( uint vertices[], Facet & facet ) const override;

  /// Return if mesh connectivities require ordering
  bool connectivity_needs_ordering( uint d0, uint d1 ) const override;

  /// Initialize mesh connectivities required by ordering
  void initialize_connectivities( Mesh & mesh ) const override;

  //--- REFINEMENT PATTERN ----------------------------------------------------

  /// Regular refinement of cell
  void
    refine_cell( Cell & cell, MeshEditor & editor, uint & current_cell ) const override;

  /// Number of vertices created by refinement pattern restricted to each
  /// entity of given topological dimensions
  uint num_refined_vertices( uint dim ) const override;

  /// Number of cells created by refinement pattern
  uint num_refined_cells() const override;

  //---------------------------------------------------------------------------

  /// Compute volume of tetrahedron
  real volume( MeshEntity const & entity ) const override;

  /// Compute diameter of tetrahedron
  real diameter( MeshEntity const & entity ) const override;

  /// Compute circumradius of tetrahedron
  real circumradius( MeshEntity const & entity ) const override;

  /// Compute inradius of interval
  real inradius( MeshEntity const & entity ) const override;

  /// Compute coordinates of midpoint
  void midpoint( MeshEntity const & entity, real * p ) const override;

  /// Compute of given facet with respect to the cell
  void normal( Cell const & cell, uint facet, real * n ) const override;

  /// Compute the area/length of given facet with respect to the cell
  real facet_area( Cell const & cell, uint facet ) const override;

  /// Check if point p intersects the entity
  bool intersects( MeshEntity const & e, Point const & p ) const override;

  /// Check if points line connecting p1 and p2 cuts the entity
  bool intersects( MeshEntity const & e,
                   Point const & p1, Point const & p2 ) const override;

  //--- REFERENCE CELL --------------------------------------------------------

  /// Create a mesh consisting of the reference cell
  void create_reference_cell( Mesh & mesh ) const override;

  /// Return coordinates of vertices in the reference cell
  real const * reference_vertex( uint i ) const override;

  //---------------------------------------------------------------------------

  /// Return description of cell type
  std::string description() const override;

  /// Display information
  void disp() const override;

  /// Check
  bool check( Cell & cell ) const override;

private:
  // Find local index of edge i according to ordering convention
  uint findEdge( uint i, Cell const & cell ) const;
};
//-----------------------------------------------------------------------------
inline uint TetrahedronCell::dim() const
{
  return 3;
}

//-----------------------------------------------------------------------------
inline uint TetrahedronCell::num_entities( uint dim ) const
{
  dolfin_assert( dim <= TD );
  return NE[3][dim];
}

//-----------------------------------------------------------------------------
inline uint TetrahedronCell::num_entities( uint d0, uint d1 ) const
{
  dolfin_assert( d0 <= TD );
  dolfin_assert( d1 <= TD );
  return NE[d0][d1];
}

//-----------------------------------------------------------------------------
inline uint TetrahedronCell::num_vertices( uint dim ) const
{
  dolfin_assert( dim <= TD );
  return NE[dim][0];
}

//-----------------------------------------------------------------------------
inline uint TetrahedronCell::orientation( Cell const & cell ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = cell.mesh().geometry();
  Array< uint > const & vertices = cell.entities( 0 );

  real const * v0 = geometry.x( vertices[0] );
  real const * v1 = geometry.x( vertices[1] );
  real const * v2 = geometry.x( vertices[2] );
  real const * v3 = geometry.x( vertices[3] );

  real a = +( ( v1[1] - v0[1] ) * ( v2[2] - v0[2] )
              - ( v1[2] - v0[2] ) * ( v2[1] - v0[1] ) )
             * ( v3[0] - v0[0] )
           + ( ( v1[2] - v0[2] ) * ( v2[0] - v0[0] )
               - ( v1[0] - v0[0] ) * ( v2[2] - v0[2] ) )
               * ( v3[1] - v0[1] )
           + ( ( v1[0] - v0[0] ) * ( v2[1] - v0[1] )
               - ( v1[1] - v0[1] ) * ( v2[0] - v0[0] ) )
               * ( v3[2] - v0[2] );

  return ( a < 0.0 ? 1 : 0 );
}

//-----------------------------------------------------------------------------
inline bool TetrahedronCell::connectivity_needs_ordering( uint d0,
                                                          uint d1 ) const
{
  dolfin_assert( d0 <= TD && d1 <= TD );
  return ( d0 > 0 && d0 > d1 );
}

//-----------------------------------------------------------------------------
inline void TetrahedronCell::initialize_connectivities( Mesh & mesh ) const
{
  mesh.init( 1, 0 );
  mesh.init( 2, 0 );
  mesh.init( 2, 1 );
  mesh.init( 3, 0 );
  mesh.init( 3, 1 );
  mesh.init( 3, 2 );
}

//-----------------------------------------------------------------------------
inline uint TetrahedronCell::num_refined_cells() const
{
  return 8;
}

//-----------------------------------------------------------------------------
inline uint TetrahedronCell::num_refined_vertices( uint dim ) const
{
  dolfin_assert( dim <= TD );
  return ( dim > 1 ? 0 : 1 );
}

//-----------------------------------------------------------------------------
inline real TetrahedronCell::volume( MeshEntity const & entity ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[3][0] );

  // Get the coordinates of the four vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );
  real const * x3 = geometry.x( vertices[3] );

  // Formula for volume from http://mathworld.wolfram.com
  real V = ( +x0[0]
               * ( x1[1] * x2[2] + x3[1] * x1[2] + x2[1] * x3[2] - x2[1] * x1[2]
                   - x1[1] * x3[2] - x3[1] * x2[2] )
             - x1[0]
                 * ( x0[1] * x2[2] + x3[1] * x0[2] + x2[1] * x3[2]
                     - x2[1] * x0[2] - x0[1] * x3[2] - x3[1] * x2[2] )
             + x2[0]
                 * ( x0[1] * x1[2] + x3[1] * x0[2] + x1[1] * x3[2]
                     - x1[1] * x0[2] - x0[1] * x3[2] - x3[1] * x1[2] )
             - x3[0]
                 * ( x0[1] * x1[2] + x1[1] * x2[2] + x2[1] * x0[2]
                     - x1[1] * x0[2] - x2[1] * x1[2] - x0[1] * x2[2] ) );

  return std::abs( V ) / 6.0;
}

//-----------------------------------------------------------------------------
inline real TetrahedronCell::diameter( MeshEntity const & entity ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[3][0] );

  // Get the coordinates of the four vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );
  real const * x3 = geometry.x( vertices[3] );

  // Compute edge lengths
  real a  = 0.0;
  real b  = 0.0;
  real c  = 0.0;
  real aa = 0.0;
  real bb = 0.0;
  real cc = 0.0;
  for ( uint i = 0; i < geometry.dim(); ++i )
  {
    a += ( x1[i] - x2[i] ) * ( x1[i] - x2[i] );
    b += ( x0[i] - x2[i] ) * ( x0[i] - x2[i] );
    c += ( x0[i] - x1[i] ) * ( x0[i] - x1[i] );
    aa += ( x0[i] - x3[i] ) * ( x0[i] - x3[i] );
    bb += ( x1[i] - x3[i] ) * ( x1[i] - x3[i] );
    cc += ( x2[i] - x3[i] ) * ( x2[i] - x3[i] );
  }

  real hmax = a;
  hmax      = std::max( b, hmax );
  hmax      = std::max( c, hmax );
  hmax      = std::max( aa, hmax );
  hmax      = std::max( bb, hmax );
  hmax      = std::max( cc, hmax );
  return std::sqrt( hmax );
}

//-----------------------------------------------------------------------------
inline real TetrahedronCell::circumradius( MeshEntity const & entity ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[3][0] );

  // Get the coordinates of the four vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );
  real const * x3 = geometry.x( vertices[3] );

  // Compute edge lengths
  real a  = 0.0;
  real b  = 0.0;
  real c  = 0.0;
  real aa = 0.0;
  real bb = 0.0;
  real cc = 0.0;
  for ( uint i = 0; i < geometry.dim(); ++i )
  {
    a += ( x1[i] - x2[i] ) * ( x1[i] - x2[i] );
    b += ( x0[i] - x2[i] ) * ( x0[i] - x2[i] );
    c += ( x0[i] - x1[i] ) * ( x0[i] - x1[i] );
    aa += ( x0[i] - x3[i] ) * ( x0[i] - x3[i] );
    bb += ( x1[i] - x3[i] ) * ( x1[i] - x3[i] );
    cc += ( x2[i] - x3[i] ) * ( x2[i] - x3[i] );
  }

  // Compute "area" of triangle with strange side lengths
  real la   = a * aa;
  real lb   = b * bb;
  real lc   = c * cc;
  real s    = 0.5 * ( la + lb + lc );
  real area = std::sqrt( s * ( s - la ) * ( s - lb ) * ( s - lc ) );

  // Formula for volume from http://mathworld.wolfram.com
  real V = ( +x0[0]
               * ( x1[1] * x2[2] + x3[1] * x1[2] + x2[1] * x3[2] - x2[1] * x1[2]
                   - x1[1] * x3[2] - x3[1] * x2[2] )
             - x1[0]
                 * ( x0[1] * x2[2] + x3[1] * x0[2] + x2[1] * x3[2]
                     - x2[1] * x0[2] - x0[1] * x3[2] - x3[1] * x2[2] )
             + x2[0]
                 * ( x0[1] * x1[2] + x3[1] * x0[2] + x1[1] * x3[2]
                     - x1[1] * x0[2] - x0[1] * x3[2] - x3[1] * x1[2] )
             - x3[0]
                 * ( x0[1] * x1[2] + x1[1] * x2[2] + x2[1] * x0[2]
                     - x1[1] * x0[2] - x2[1] * x1[2] - x0[1] * x2[2] ) );

  // Formula for circumradius from http://mathworld.wolfram.com
  return area / ( 6.0 * V );
}

//-----------------------------------------------------------------------------
inline real TetrahedronCell::inradius( MeshEntity const & entity ) const
{
  return circumradius( entity ) / 3.0;
}

//-----------------------------------------------------------------------------
inline void TetrahedronCell::midpoint( MeshEntity const & entity,
                                       real *             p ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[3][0] );

  // Get the coordinates of the vertices
  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );
  real const * x3 = geometry.x( vertices[3] );

  for ( uint d = 0; d < geometry.dim(); ++d )
  {
    p[d] = 0.25 * ( x0[d] + x1[d] + x2[d] + x3[d] );
  }
}

//-----------------------------------------------------------------------------
inline void
  TetrahedronCell::normal( Cell const & cell, uint facet, real * n ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  // Create facet from the mesh and local facet number
  Cell &               c = const_cast< Cell & >( cell );
  Facet                f( c.mesh(), c.entities( 2 )[facet] );
  MeshGeometry const & geometry = cell.mesh().geometry();

  // Get coordinates of opposite vertex
  real const * p0 = geometry.x( cell.entities( 0 )[facet] );

  // Get coordinates of facet vertices
  Array< uint > const & vertices = f.entities( 0 );

  real const * p1 = geometry.x( vertices[0] );
  real const * p2 = geometry.x( vertices[1] );
  real const * p3 = geometry.x( vertices[2] );

  // n = e1 ^ e2
  n[0] = ( p2[1] - p1[1] ) * ( p3[2] - p1[2] )
         - ( p2[2] - p1[2] ) * ( p3[1] - p1[1] );
  n[1] = ( p2[2] - p1[2] ) * ( p3[0] - p1[0] )
         - ( p2[0] - p1[0] ) * ( p3[2] - p1[2] );
  n[2] = ( p2[0] - p1[0] ) * ( p3[1] - p1[1] )
         - ( p2[1] - p1[1] ) * ( p3[0] - p1[0] );
  real const nn = std::sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2] );
  n[0] /= nn;
  n[1] /= nn;
  n[2] /= nn;
  //
  if ( n[0] * ( p1[0] - p0[0] ) + n[1] * ( p1[1] - p0[1] )
         + n[2] * ( p1[2] - p0[2] )
       < 0.0 )
  {
    n[0] *= -1.0;
    n[1] *= -1.0;
    n[2] *= -1.0;
  }
}

//-----------------------------------------------------------------------------
inline real TetrahedronCell::facet_area( Cell const & cell, uint facet ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  // Create facet from the mesh and local facet number
  Cell & c = const_cast< Cell & >( cell );
  Facet  f( c.mesh(), c.entities( 2 )[facet] );

  // Get the coordinates of the three vertices
  MeshGeometry const &  geometry = cell.mesh().geometry();
  Array< uint > const & vertices = f.entities( 0 );

  real const * x0 = geometry.x( vertices[0] );
  real const * x1 = geometry.x( vertices[1] );
  real const * x2 = geometry.x( vertices[2] );

  // Compute area of triangle embedded in R^3
  real v0 = ( x0[1] * x1[2] + x0[2] * x2[1] + x1[1] * x2[2] )
            - ( x2[1] * x1[2] + x2[2] * x0[1] + x1[1] * x0[2] );
  real v1 = ( x0[2] * x1[0] + x0[0] * x2[2] + x1[2] * x2[0] )
            - ( x2[2] * x1[0] + x2[0] * x0[2] + x1[2] * x0[0] );
  real v2 = ( x0[0] * x1[1] + x0[1] * x2[0] + x1[0] * x2[1] )
            - ( x2[0] * x1[1] + x2[1] * x0[0] + x1[0] * x0[1] );

  // Formula for area from http://mathworld.wolfram.com
  return 0.5 * std::sqrt( v0 * v0 + v1 * v1 + v2 * v2 );
}

//-----------------------------------------------------------------------------
inline real const * TetrahedronCell::reference_vertex( uint i ) const
{
  return &VC[i][0];
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_TETRAHEDRON_CELL_H */
