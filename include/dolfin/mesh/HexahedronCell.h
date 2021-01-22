// Copyright (C) 2015 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_HEXAHEDRON_CELL_H
#define __DOLFIN_HEXAHEDRON_CELL_H

#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Edge.h>
#include <dolfin/mesh/entities/Facet.h>

namespace dolfin
{

/**
 *  @class  HexahedronCell
 *
 *  @brief  This class implements functionality for hexahedral meshes.
 *
 */

class HexahedronCell : public CellType
{
  // UFC: Topological Dimension
  static size_t const TD = 3;

  // UFC: Number of Entities
  static size_t const NE[4][4];

  // UFC: Vertex Coordinates
  static real const VC[8][3];

  // UFC: Edge - Incident Vertices
  static size_t const EIV[12][2];

  // UFC: Edge - Non-Incident Vertices
  static size_t const ENV[12][6];

  // UFC: Face - Incident Vertices
  static size_t const FIV[6][4];

  // UFC: Face - Non-Incident Vertices
  static size_t const FNV[6][4];

public:
  /// Specify cell type and facet type
  HexahedronCell();

  ///
  ~HexahedronCell() override = default;

  /// Clone pattern
  auto clone() const -> CellType * override
  {
    return new HexahedronCell( *this );
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

  /// Order entities locally (connectivity 1-0, 2-0, 2-1, 3-0, 3-1, 3-2)
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

  /// Compute (generalized) volume (area) of hexahedron
  auto volume( MeshEntity const & entity ) const -> real override;

  /// Compute diameter of hexahedron
  auto diameter( MeshEntity const & entity ) const -> real override;

  /// Compute circumradius of hexahedron
  auto circumradius( MeshEntity const & entity ) const -> real override;

  /// Compute inradius of hexahedron
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

  // Find local index of face i according to ordering convention
  auto findFace( size_t i, Cell const & cell ) const -> size_t;
};

//-----------------------------------------------------------------------------
inline auto HexahedronCell::dim() const -> size_t
{
  return 3;
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::num_entities( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= TD );
  return NE[3][dim];
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::num_entities( size_t d0, size_t d1 ) const -> size_t
{
  dolfin_assert( d0 <= TD );
  dolfin_assert( d1 <= TD );
  return NE[d0][d1];
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::num_vertices( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= TD );
  return NE[dim][0];
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::orientation( Cell const & cell ) const -> size_t
{
  dolfin_assert( cell.type() == this->cell_type );

  // Get the coordinates of vertices v0, v1, v3 and v4
  MeshGeometry const &    geometry = cell.mesh().geometry();
  Array< size_t > const & vertices = cell.entities( 0 );
  real const *            v0       = geometry.x( vertices[0] );
  real const *            v1       = geometry.x( vertices[1] );
  real const *            v3       = geometry.x( vertices[3] );
  real const *            v4       = geometry.x( vertices[4] );

  // Check whether (v0v1, v0v3, v0v4) is counter-clockwise
  real const S = +( ( v1[1] - v0[1] ) * ( v3[2] - v0[2] )
                    - ( v1[2] - v0[2] ) * ( v3[1] - v0[1] ) )
                   * ( v4[0] - v0[0] )
                 + ( ( v1[2] - v0[2] ) * ( v3[0] - v0[0] )
                     - ( v1[0] - v0[0] ) * ( v3[2] - v0[2] ) )
                     * ( v4[1] - v0[1] )
                 + ( ( v1[0] - v0[0] ) * ( v3[1] - v0[1] )
                     - ( v1[1] - v0[1] ) * ( v3[0] - v0[0] ) )
                     * ( v4[2] - v0[2] );
  return ( S < 0.0 ? 1 : 0 );
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::connectivity_needs_ordering( size_t d0,
                                                         size_t d1 ) const
  -> bool
{
  dolfin_assert( d0 <= TD && d1 <= TD );
  // Do not order cell - vertices connectivities
  return ( d0 > 0 && d0 > d1 ) && !( d0 == TD && d1 == 0 );
}

//-----------------------------------------------------------------------------
inline void HexahedronCell::initialize_connectivities( Mesh & mesh ) const
{
  mesh.init( 1, 0 );
  mesh.init( 2, 0 );
  mesh.init( 2, 1 );
  mesh.init( 3, 0 );
  mesh.init( 3, 1 );
  mesh.init( 3, 2 );
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::num_refined_cells() const -> size_t
{
  return 8;
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::num_refined_vertices( size_t dim ) const -> size_t
{
  dolfin_assert( dim <= TD );
  MAYBE_UNUSED( dim );
  return 1;
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::volume( MeshEntity const & entity ) const -> real
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[3][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &    geometry = entity.mesh().geometry();
  Array< size_t > const & vertices = entity.entities( 0 );
  real const *            v0       = geometry.x( vertices[0] );
  real const *            v1       = geometry.x( vertices[1] );
  real const *            v3       = geometry.x( vertices[3] );
  real const *            v4       = geometry.x( vertices[4] );

  // Compute mixed product of (v0v1, v0v3, v0v4)
  real const S = +( ( v1[1] - v0[1] ) * ( v3[2] - v0[2] )
                    - ( v1[2] - v0[2] ) * ( v3[1] - v0[1] ) )
                   * ( v4[0] - v0[0] )
                 + ( ( v1[2] - v0[2] ) * ( v3[0] - v0[0] )
                     - ( v1[0] - v0[0] ) * ( v3[2] - v0[2] ) )
                     * ( v4[1] - v0[1] )
                 + ( ( v1[0] - v0[0] ) * ( v3[1] - v0[1] )
                     - ( v1[1] - v0[1] ) * ( v3[0] - v0[0] ) )
                     * ( v4[2] - v0[2] );
  return std::fabs( S );
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::diameter( MeshEntity const & entity ) const -> real
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[3][0] );

  // Get the coordinates of the three vertices
  MeshGeometry const &    geometry = entity.mesh().geometry();
  Array< size_t > const & vertices = entity.entities( 0 );
  real const *            x0       = geometry.x( vertices[0] );
  real const *            x1       = geometry.x( vertices[1] );
  real const *            x2       = geometry.x( vertices[2] );
  real const *            x3       = geometry.x( vertices[3] );
  real const *            x4       = geometry.x( vertices[4] );
  real const *            x5       = geometry.x( vertices[5] );
  real const *            x6       = geometry.x( vertices[6] );
  real const *            x7       = geometry.x( vertices[7] );

  // Compute maximum diagonal
  real d0 = 0.0;
  real d1 = 0.0;
  real d2 = 0.0;
  real d3 = 0.0;
  for ( size_t i = 0; i < geometry.dim(); ++i )
  {
    d0 += ( x6[i] - x0[i] ) * ( x6[i] - x0[i] );
    d1 += ( x7[i] - x1[i] ) * ( x7[i] - x1[i] );
    d2 += ( x4[i] - x2[i] ) * ( x4[i] - x2[i] );
    d3 += ( x5[i] - x3[i] ) * ( x5[i] - x3[i] );
  }
  return std::sqrt( std::max( std::max( d0, d1 ), std::max( d2, d3 ) ) );
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::circumradius( MeshEntity const & entity ) const
  -> real
{
  /// @todo No better idea for now
  return 0.5 * std::sqrt( 3.0 ) * this->diameter( entity );
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::inradius( MeshEntity const & entity ) const -> real
{
  /// @todo No better idea for now
  return 0.5 * this->diameter( entity );
}

//-----------------------------------------------------------------------------
inline void HexahedronCell::midpoint( MeshEntity const & entity,
                                      real *             p ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[3][0] );

  // Get the coordinates of the vertices
  MeshGeometry const &    geometry = entity.mesh().geometry();
  Array< size_t > const & vertices = entity.entities( 0 );
  real const *            x0       = geometry.x( vertices[0] );
  real const *            x1       = geometry.x( vertices[1] );
  real const *            x2       = geometry.x( vertices[2] );
  real const *            x3       = geometry.x( vertices[3] );
  real const *            x4       = geometry.x( vertices[4] );
  real const *            x5       = geometry.x( vertices[5] );
  real const *            x6       = geometry.x( vertices[6] );
  real const *            x7       = geometry.x( vertices[7] );
  for ( size_t i = 0; i < geometry.dim(); ++i )
  {
    p[i] =
      0.125 * ( x0[i] + x1[i] + x2[i] + x3[i] + x4[i] + x5[i] + x6[i] + x7[i] );
  }
}

//-----------------------------------------------------------------------------
inline void
  HexahedronCell::normal( Cell const & cell, size_t facet, real * n ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  // Create facet from the mesh and local facet number
  Cell &               c = const_cast< Cell & >( cell );
  Facet                f( c.mesh(), c.entities( 1 )[facet] );
  MeshGeometry const & geometry = cell.mesh().geometry();

  // first non-incident vertex
  real const * vf = geometry.x( c.entities( 0 )[FNV[facet][0]] );

  // vertices on the facet
  Array< size_t > const & vertices = f.entities( 0 );
  real const *            v0       = geometry.x( vertices[FIV[facet][0]] );
  real const *            v1       = geometry.x( vertices[FIV[facet][1]] );
  real const *            v3       = geometry.x( vertices[FIV[facet][3]] );

  // Vector normal to facet
  n[0] = +( v1[1] - v0[1] ) * ( v3[2] - v0[2] )
         - ( v1[2] - v0[2] ) * ( v3[1] - v0[1] );
  n[1] = +( v1[2] - v0[2] ) * ( v3[0] - v0[0] )
         - ( v1[0] - v0[0] ) * ( v3[2] - v0[2] );
  n[2] = +( v1[0] - v0[0] ) * ( v3[1] - v0[1] )
         - ( v1[1] - v0[1] ) * ( v3[0] - v0[0] );
  real const nn = std::sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2] );
  n[0] /= nn;
  n[1] /= nn;
  n[2] /= nn;

  // Flip direction of normal so it points outward
  if ( +n[0] * ( vf[0] - v0[0] ) + n[1] * ( vf[1] - v0[1] )
         + n[2] * ( vf[2] - v0[2] )
       > 0.0 )
  {
    n[0] *= -1.0;
    n[1] *= -1.0;
    n[2] *= -1.0;
  }
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::facet_area( Cell const & cell, size_t facet ) const
  -> real
{
  dolfin_assert( cell.type() == this->cell_type );

  // Create facet from the mesh and local facet number
  Cell & c = const_cast< Cell & >( cell );
  Facet  f( c.mesh(), c.entities( 1 )[facet] );

  // Get the coordinates of the two vertices
  MeshGeometry const & geometry = cell.mesh().geometry();
  real const *         p        = geometry.x( f.entities( 0 )[FIV[facet][0]] );
  real const *         q        = geometry.x( f.entities( 0 )[FIV[facet][1]] );
  real const *         r        = geometry.x( f.entities( 0 )[FIV[facet][2]] );
  real const *         s        = geometry.x( f.entities( 0 )[FIV[facet][3]] );

  // Compute the area with diagonal formula
  real c0 =
    ( r[1] - p[1] ) * ( s[2] - q[2] ) - ( r[2] - p[2] ) * ( s[1] - q[1] );
  real c1 =
    ( r[2] - p[2] ) * ( s[0] - q[0] ) - ( r[0] - p[0] ) * ( s[2] - q[2] );
  real c2 =
    ( r[0] - p[0] ) * ( s[1] - q[1] ) - ( r[1] - p[1] ) * ( s[0] - q[0] );
  return 0.5 * std::sqrt( c0 * c0 + c1 * c1 + c2 * c2 );
}

//-----------------------------------------------------------------------------
inline auto HexahedronCell::reference_vertex( size_t i ) const -> real const *
{
  return &VC[i][0];
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_HEXAHEDRON_CELL_H */
