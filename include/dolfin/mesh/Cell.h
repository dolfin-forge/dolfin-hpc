// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CELL_H
#define __DOLFIN_CELL_H

#include "Point.h"
#include "CellType.h"
#include "MeshEntity.h"
#include "MeshEntityIterator.h"

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
  Cell(Mesh& mesh, uint index) :
      MeshEntity(mesh, mesh.topology_dimension(), index)
  {
  }

  /// Destructor
  ~Cell()
  {
  }

  /// Return type of cell
  inline CellType::Type type() const
  {
    return mesh_.type().cellType();
  }

  /// Compute orientation of cell (0 is right, 1 is left)
  inline real orientation() const
  {
    return mesh_.type().orientation(*this);
  }

  /// Compute (generalized) volume of cell
  inline real volume() const
  {
    return mesh_.type().volume(*this);
  }

  /// Compute diameter of cell
  inline real diameter() const
  {
    return mesh_.type().diameter(*this);
  }

  /// Compute circumradius of cell
  inline real circumradius() const
  {
    return mesh_.type().circumradius(*this);
  }

  /// Compute circumradius of cell
  inline real inradius() const
  {
    return mesh_.type().inradius(*this);
  }

  /// Compute normal of given facet with respect to the cell
  inline Point normal(uint facet) const
  {
    Point n;
    mesh_.type().normal(*this, facet, &n[0]);
    return n;
  }

  /// Compute normal of given facet with respect to the cell
  inline void normal(uint facet, real * n) const
  {
    mesh_.type().normal(*this, facet, n);
  }

  /// Compute the area/length of given facet with respect to the cell
  inline real facet_area(uint facet) const
  {
    return mesh_.type().facet_area(*this, facet);
  }

  /// Compute coordinates of cell midpoint
  inline Point midpoint() const
  {
    Point p;
    mesh_.type().midpoint(*this, &p[0]);
    return p;
  }

  /// Compute coordinates of cell midpoint
  inline void midpoint(real * p) const
  {
    mesh_.type().midpoint(*this, p);
  }

  //--- ITERATOR --------------------------------------------------------------

  typedef CellIterator iterator;

  struct shared : SharedIterator
  {
    shared(Mesh& M) : SharedIterator(M.topology().distdata()[M.type().dim()]) {}
    shared(MeshTopology& T) : SharedIterator(T.distdata()[T.dim()]) {}
  };

  struct ghost : GhostIterator
  {
    ghost(Mesh& M) : GhostIterator(M.topology().distdata()[M.type().dim()]) {}
    ghost(MeshTopology& T) : GhostIterator(T.distdata()[T.dim()]) {}
  };

  struct owned : OwnedIterator
  {
    owned(Mesh& M) : OwnedIterator(M.topology().distdata()[M.type().dim()]) {}
    owned(MeshTopology& T) : OwnedIterator(T.distdata()[T.dim()]) {}
  };

  //--- Entity relation -------------------------------------------------------

  typedef Face lower_dimensional;

};

/**
 *  @class  CellIterator
 *
 *  @brief  A CellIterator is a MeshEntityIterator of topological codimension 0.
 *
 */

class CellIterator : public MeshEntityIterator
{
public:

  CellIterator(Mesh& mesh) :
      MeshEntityIterator(mesh, mesh.topology_dimension())
  {
  }

  CellIterator(MeshEntity& entity) :
      MeshEntityIterator(entity, entity.mesh().topology_dimension())
  {
  }

  inline Cell* operator->()
  {
    return static_cast<Cell*>(MeshEntityIterator::operator->());
  }

  inline Cell& operator*()
  {
    return *operator->();
  }

  inline Cell& operator[](uint i)
  {
    return static_cast<Cell&>(MeshEntityIterator::operator[](i));
  }

};

} /* namespace dolfin */

#endif /* __DOLFIN_CELL_H */
