// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_INTERVAL_CELL_H
#define __DOLFIN_INTERVAL_CELL_H

#include <dolfin/common/constants.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

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
  static uint const TD = 1;

  // UFC: Number of Entities
  static uint const NE[2][2];

  // UFC: Vertex Coordinates
  static real const VC[2][1];

  // UFC: Edge - Incident Vertices
  static uint const EIV[1][2];

public:
  /// Specify cell type and facet type
  IntervalCell();

  ///
  ~IntervalCell() override = default;

  /// Clone pattern
  CellType * clone() const override
  {
    return new IntervalCell( *this );
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

  /// Order entities locally (connectivity 1-0)
  void order_entities( MeshTopology & topology, uint i ) const override;

  /// Order vertices such that the facet is right-oriented w.r.t. facet normal
  void order_facet( uint vertices[], Facet & facet ) const override;

  /// Return if mesh connectivities require ordering
  bool connectivity_needs_ordering( uint d0, uint d1 ) const override;

  /// Initialize mesh connectivities required by ordering
  void initialize_connectivities( Mesh & mesh ) const override;

  //--- REFINEMENT PATTERN ----------------------------------------------------

  /// Refine cell uniformly
  void
    refine_cell( Cell & cell, MeshEditor & editor, uint & current_cell ) const override;

  /// Number of cells created by refinement pattern
  uint num_refined_cells() const override;

  /// Number of vertices created by refinement pattern restricted to each
  /// entity of given topological dimensions
  uint num_refined_vertices( uint dim ) const override;

  //---------------------------------------------------------------------------

  /// Compute (generalized) volume (length) of interval
  real volume( MeshEntity const & entity ) const override;

  /// Compute diameter of interval
  real diameter( MeshEntity const & entity ) const override;

  /// Compute circumradius of interval
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
                   Point const &      p1,
                   Point const &      p2 ) const override;

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
};

//-----------------------------------------------------------------------------
inline uint IntervalCell::dim() const
{
  return 1;
}

//-----------------------------------------------------------------------------
inline uint IntervalCell::num_entities( uint dim ) const
{
  dolfin_assert( dim <= TD );
  return NE[1][dim];
}

//-----------------------------------------------------------------------------
inline uint IntervalCell::num_entities( uint d0, uint d1 ) const
{
  dolfin_assert( d0 <= TD );
  dolfin_assert( d1 <= TD );
  return NE[d0][d1];
}

//-----------------------------------------------------------------------------
inline uint IntervalCell::num_vertices( uint dim ) const
{
  dolfin_assert( dim <= TD );
  return NE[dim][0];
}

//-----------------------------------------------------------------------------
inline uint IntervalCell::orientation( Cell const & cell ) const
{
  dolfin_assert( cell.type() == this->cell_type );
  Point v01 =   Point( cell.entities( 0 )[1], 0.0, 0.0 )
              - Point( cell.entities( 0 )[0], 0.0, 0.0 );
  Point n( -v01[1], v01[0], 0.0 );

  return ( n.dot( v01 ) < 0.0 ? 1 : 0 );
}

//-----------------------------------------------------------------------------
inline void IntervalCell::order_facet( uint[], Facet & ) const
{
  // Do nothing
}

//-----------------------------------------------------------------------------
inline bool IntervalCell::connectivity_needs_ordering( uint d0, uint d1 ) const
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
inline uint IntervalCell::num_refined_cells() const
{
  return 2;
}

//-----------------------------------------------------------------------------
inline uint IntervalCell::num_refined_vertices( uint dim ) const
{
  dolfin_assert( dim <= TD );
  MAYBE_UNUSED( dim );
  return 1;
}

//-----------------------------------------------------------------------------
inline real IntervalCell::volume( MeshEntity const & entity ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[1][0] );

  // Get mesh geometry
  MeshGeometry const & geometry = entity.mesh().geometry();

  // Get the coordinates of the two vertices
  Array< uint > const & vertices = entity.entities( 0 );
  real const *          x0       = geometry.x( vertices[0] );
  real const *          x1       = geometry.x( vertices[1] );

  // Compute length of interval (line segment)
  real sum = 0.0;
  for ( uint i = 0; i < geometry.dim(); ++i )
  {
    sum += ( x1[i] - x0[i] ) * ( x1[i] - x0[i] );
  }

  return std::sqrt( sum );
}

//-----------------------------------------------------------------------------
inline real IntervalCell::diameter( MeshEntity const & entity ) const
{
  // Diameter is same as volume for interval (line segment)
  return volume( entity );
}

//-----------------------------------------------------------------------------
inline real IntervalCell::circumradius( MeshEntity const & entity ) const
{
  // Circumradius is same as volume for interval (line segment)
  return volume( entity );
}

//-----------------------------------------------------------------------------
inline real IntervalCell::inradius( MeshEntity const & entity ) const
{
  // Inradius is same as volume for interval (line segment)
  return volume( entity );
}

//-----------------------------------------------------------------------------
inline void IntervalCell::midpoint( MeshEntity const & entity, real * p ) const
{
  dolfin_assert( entity.dim() == TD );
  dolfin_assert( entity.num_entities( 0 ) == NE[1][0] );

  MeshGeometry const &  geometry = entity.mesh().geometry();
  Array< uint > const & vertices = entity.entities( 0 );
  real const *          x0       = geometry.x( vertices[0] );
  real const *          x1       = geometry.x( vertices[1] );
  uint const            gdim     = geometry.dim();
  for ( uint d = 0; d < gdim; ++d )
  {
    p[d] = 0.5 * ( x0[d] + x1[d] );
  }
}

//-----------------------------------------------------------------------------
inline void
  IntervalCell::normal( Cell const & cell, uint facet, real * n ) const
{
  dolfin_assert( cell.type() == this->cell_type );

  MeshGeometry const &  geometry = cell.mesh().geometry();
  Array< uint > const & vertices = cell.entities( 0 );

  Point p0 = geometry.point( vertices[facet] );
  Point p1 = geometry.point( vertices[( facet + 1 ) % 2] );
  real nn  = p0.dist( p1 );

  for ( uint d = 0; d < geometry.dim(); ++d )
  {
    n[d] /= nn;
  }
}

//-----------------------------------------------------------------------------
inline bool IntervalCell::intersects( MeshEntity const & e,
                                      Point const &      p ) const
{
  dolfin_assert( e.dim() == TD );
  dolfin_assert( e.num_entities( 0 ) == NE[1][0] );

  // Get the coordinates of the vertices
  MeshGeometry const &  geometry = e.mesh().geometry();
  Array< uint > const & vertices = e.entities( 0 );

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
inline real IntervalCell::facet_area( Cell const & cell, uint ) const
{
  dolfin_assert( cell.type() == this->cell_type );
  MAYBE_UNUSED( cell );
  return 0.0;
}

//-----------------------------------------------------------------------------
inline real const * IntervalCell::reference_vertex( uint i ) const
{
  return &VC[i][0];
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_INTERVAL_CELL_H */
