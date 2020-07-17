// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CELL_H
#define __DOLFIN_CELL_H

#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/Point.h>

#include <dolfin/common/GhostIterator.h>
#include <dolfin/common/OwnedIterator.h>
#include <dolfin/common/SharedIterator.h>

namespace dolfin
{

/**
 *  @class  Cell
 *
 *  @brief  A Cell is a MeshEntity of topological codimension 0.
 *
 */

class CellIterator;

class Cell : public MeshEntity
{

public:
  /// Constructor
  Cell( Mesh & mesh, uint index )
    : MeshEntity( mesh, mesh.topology_dimension(), index )
  {
  }

  /// Destructor
  ~Cell()
  {
  }

  /// Return type of cell
  inline CellType::Type type() const;

  /// Compute orientation of cell (0 is right, 1 is left)
  inline real orientation() const;

  /// Compute (generalized) volume of cell
  inline real volume() const;

  /// Compute diameter of cell
  inline real diameter() const;

  /// Compute circumradius of cell
  inline real circumradius() const;

  /// Compute circumradius of cell
  inline real inradius() const;

  /// Compute normal of given facet with respect to the cell
  inline Point normal( uint facet ) const;

  /// Compute normal of given facet with respect to the cell
  inline void normal( uint facet, real * n ) const;

  /// Compute the area/length of given facet with respect to the cell
  inline real facet_area( uint facet ) const;

  /// Compute coordinates of cell midpoint
  inline Point midpoint() const;

  /// Compute coordinates of cell midpoint
  inline void midpoint( real * p ) const;

  //--- ITERATOR --------------------------------------------------------------

  using iterator = CellIterator;

  struct shared : SharedIterator
  {
    shared( Mesh & M )
      : SharedIterator( M.topology().distdata()[M.type().dim()] )
    {
    }
    shared( MeshTopology & T )
      : SharedIterator( T.distdata()[T.dim()] )
    {
    }
  };

  struct ghost : GhostIterator
  {
    ghost( Mesh & M )
      : GhostIterator( M.topology().distdata()[M.type().dim()] )
    {
    }
    ghost( MeshTopology & T )
      : GhostIterator( T.distdata()[T.dim()] )
    {
    }
  };

  struct owned : OwnedIterator
  {
    owned( Mesh & M )
      : OwnedIterator( M.topology().distdata()[M.type().dim()] )
    {
    }
    owned( MeshTopology & T )
      : OwnedIterator( T.distdata()[T.dim()] )
    {
    }
  };

  //--- Entity relation -------------------------------------------------------

  using lower_dimensional = Face;
};

//-----------------------------------------------------------------------------
inline CellType::Type Cell::type() const
{
  return mesh_.type().cellType();
}

//-----------------------------------------------------------------------------
inline real Cell::orientation() const
{
  return mesh_.type().orientation( *this );
}

//-----------------------------------------------------------------------------
inline real Cell::volume() const
{
  return mesh_.type().volume( *this );
}

//-----------------------------------------------------------------------------
inline real Cell::diameter() const
{
  return mesh_.type().diameter( *this );
}

//-----------------------------------------------------------------------------
inline real Cell::circumradius() const
{
  return mesh_.type().circumradius( *this );
}

//-----------------------------------------------------------------------------
inline real Cell::inradius() const
{
  return mesh_.type().inradius( *this );
}

//-----------------------------------------------------------------------------
inline Point Cell::normal( uint facet ) const
{
  Point n;
  mesh_.type().normal( *this, facet, &n[0] );
  return n;
}

//-----------------------------------------------------------------------------
inline void Cell::normal( uint facet, real * n ) const
{
  mesh_.type().normal( *this, facet, n );
}

//-----------------------------------------------------------------------------
inline real Cell::facet_area( uint facet ) const
{
  return mesh_.type().facet_area( *this, facet );
}

//-----------------------------------------------------------------------------
inline Point Cell::midpoint() const
{
  Point p;
  mesh_.type().midpoint( *this, &p[0] );
  return p;
}

//-----------------------------------------------------------------------------
inline void Cell::midpoint( real * p ) const
{
  mesh_.type().midpoint( *this, p );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_CELL_H */
