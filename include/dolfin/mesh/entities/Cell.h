// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CELL_H
#define __DOLFIN_CELL_H

#include <dolfin/common/GhostIterator.h>
#include <dolfin/common/OwnedIterator.h>
#include <dolfin/common/SharedIterator.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/celltypes/CellType.h>
#include <dolfin/mesh/entities/MeshEntity.h>

namespace dolfin
{

class CellIterator;

//-----------------------------------------------------------------------------

/**
 *  @class  Cell
 *
 *  @brief  A Cell is a MeshEntity of topological codimension 0.
 *
 */

class Cell : public MeshEntity
{

public:
  /// Constructor
  Cell( Mesh & mesh, size_t index );

  /// Destructor
  ~Cell() = default;

  /// Return type of cell
  inline auto type() const -> CellType::Type;

  /// Compute orientation of cell (0 is right, 1 is left)
  inline auto orientation() const -> real;

  /// Compute (generalized) volume of cell
  inline auto volume() const -> real;

  /// Compute diameter of cell
  inline auto diameter() const -> real;

  /// Compute circumradius of cell
  inline auto circumradius() const -> real;

  /// Compute circumradius of cell
  inline auto inradius() const -> real;

  /// Compute normal of given facet with respect to the cell
  inline auto normal( size_t facet ) const -> Point;

  /// Compute normal of given facet with respect to the cell
  inline auto normal( size_t facet, real * n ) const -> void;

  /// Compute the area/length of given facet with respect to the cell
  inline auto facet_area( size_t facet ) const -> real;

  /// Compute coordinates of cell midpoint
  inline auto midpoint() const -> Point;

  /// Compute coordinates of cell midpoint
  inline auto midpoint( real * p ) const -> void;

  //--- ITERATOR --------------------------------------------------------------

  using iterator = CellIterator;

  struct shared : SharedIterator
  {
    shared( Mesh & M );
    shared( MeshTopology & T );
  };

  struct ghost : GhostIterator
  {
    ghost( Mesh & M );
    ghost( MeshTopology & T );
  };

  struct owned : OwnedIterator
  {
    owned( Mesh & M );
    owned( MeshTopology & T );
  };

  //--- Entity relation -------------------------------------------------------

  using lower_dimensional = Face;
};

//-----------------------------------------------------------------------------

inline Cell::Cell( Mesh & mesh, size_t index )
  : MeshEntity( mesh, mesh.topology_dimension(), index )
{
}

//-----------------------------------------------------------------------------

inline auto Cell::type() const -> CellType::Type
{
  return mesh_.type().cellType();
}

//-----------------------------------------------------------------------------

inline auto Cell::orientation() const -> real
{
  return mesh_.type().orientation( *this );
}

//-----------------------------------------------------------------------------

inline auto Cell::volume() const -> real
{
  return mesh_.type().volume( *this );
}

//-----------------------------------------------------------------------------

inline auto Cell::diameter() const -> real
{
  return mesh_.type().diameter( *this );
}

//-----------------------------------------------------------------------------

inline auto Cell::circumradius() const -> real
{
  return mesh_.type().circumradius( *this );
}

//-----------------------------------------------------------------------------

inline auto Cell::inradius() const -> real
{
  return mesh_.type().inradius( *this );
}

//-----------------------------------------------------------------------------

inline auto Cell::normal( size_t facet ) const -> Point
{
  Point n;
  mesh_.type().normal( *this, facet, &n[0] );
  return n;
}

//-----------------------------------------------------------------------------

inline auto Cell::normal( size_t facet, real * n ) const -> void
{
  mesh_.type().normal( *this, facet, n );
}

//-----------------------------------------------------------------------------

inline auto Cell::facet_area( size_t facet ) const -> real
{
  return mesh_.type().facet_area( *this, facet );
}

//-----------------------------------------------------------------------------

inline auto Cell::midpoint() const -> Point
{
  Point p;
  mesh_.type().midpoint( *this, &p[0] );
  return p;
}

//-----------------------------------------------------------------------------

inline auto Cell::midpoint( real * p ) const -> void
{
  mesh_.type().midpoint( *this, p );
}

//-----------------------------------------------------------------------------

inline Cell::shared::shared( Mesh & M )
  : SharedIterator( M.topology().distdata()[M.type().dim()] )
{
}

//-----------------------------------------------------------------------------

inline Cell::shared::shared( MeshTopology & T )
  : SharedIterator( T.distdata()[T.dim()] )
{
}

//-----------------------------------------------------------------------------

inline Cell::ghost::ghost( Mesh & M )
  : GhostIterator( M.topology().distdata()[M.type().dim()] )
{
}

//-----------------------------------------------------------------------------

inline Cell::ghost::ghost( MeshTopology & T )
  : GhostIterator( T.distdata()[T.dim()] )
{
}

//-----------------------------------------------------------------------------

inline Cell::owned::owned( Mesh & M )
  : OwnedIterator( M.topology().distdata()[M.type().dim()] )
{
}

//-----------------------------------------------------------------------------

inline Cell::owned::owned( MeshTopology & T )
  : OwnedIterator( T.distdata()[T.dim()] )
{
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_CELL_H */
