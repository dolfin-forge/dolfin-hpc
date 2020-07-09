// Copyright (C) 2007-2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_POINT_CELL_H
#define __DOLFIN_POINT_CELL_H

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
  ~PointCell();

  /// Clone pattern
  CellType * clone() const
  {
    return new PointCell( *this );
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
  void refine_cell( Cell & cell, MeshEditor & editor,
                    uint & current_cell ) const;

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

  /// Check if point p intersects the cell
  bool intersects( MeshEntity const & e, Point const & p ) const;

  /// Check if points line connecting p1 and p2 cuts the cell
  bool intersects( MeshEntity const & e,
                   Point const & p1, Point const & p2 ) const;

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
};

//-----------------------------------------------------------------------------
inline uint PointCell::dim() const
{
  return 0;
}

//-----------------------------------------------------------------------------
inline uint PointCell::num_entities( uint dim ) const
{
  dolfin_assert( dim <= TD );
  MAYBE_UNUSED( dim );
  return 1;
}

//-----------------------------------------------------------------------------
inline uint PointCell::num_entities( uint d0, uint d1 ) const
{
  dolfin_assert( d0 <= TD );
  dolfin_assert( d1 <= TD );
  MAYBE_UNUSED( d0 );
  MAYBE_UNUSED( d1 );
  return 1;
}

//-----------------------------------------------------------------------------
inline uint PointCell::num_vertices( uint dim ) const
{
  dolfin_assert( dim <= TD );
  MAYBE_UNUSED( dim );
  return 1;
}

//-----------------------------------------------------------------------------
inline uint PointCell::orientation( Cell const & ) const
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
inline bool PointCell::connectivity_needs_ordering( uint d0, uint d1 ) const
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
inline uint PointCell::num_refined_cells() const
{
  return 1;
}

//-----------------------------------------------------------------------------
inline uint PointCell::num_refined_vertices( uint ) const
{
  return 1;
}

//-----------------------------------------------------------------------------
inline real PointCell::volume( MeshEntity const & ) const
{
  return 0.0;
}

//-----------------------------------------------------------------------------
inline real PointCell::diameter( MeshEntity const & ) const
{
  return 0.0;
}

//-----------------------------------------------------------------------------
inline real PointCell::circumradius( MeshEntity const & ) const
{
  return 0.0;
}

//-----------------------------------------------------------------------------
inline real PointCell::inradius( MeshEntity const & ) const
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
inline real PointCell::facet_area( Cell const &, uint ) const
{
  return 0.0;
}

//-----------------------------------------------------------------------------
inline bool PointCell::intersects( MeshEntity const & e, Point const & p ) const
{
  return abscmp( p.dist( e.mesh().geometry().point( e.index() ) ), 0.0 );
}

//-----------------------------------------------------------------------------
inline bool PointCell::intersects( MeshEntity const &,
                                   Point const &,
                                   Point const & ) const
{
  error( "PointCell::intersects() not implemented." );
  return true;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_POINT_CELL_H */
