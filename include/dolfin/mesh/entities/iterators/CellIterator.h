// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CELL_ITERATOR_H
#define __DOLFIN_CELL_ITERATOR_H

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/MeshEntity.h>
#include <dolfin/mesh/entities/iterators/MeshEntityIterator.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

/**
 *  @class  CellIterator
 *
 *  @brief  A CellIterator is a MeshEntityIterator of topological codimension 0.
 *
 */

class CellIterator : public MeshEntityIterator
{
public:
  CellIterator( Mesh & mesh );
  CellIterator( MeshEntity & entity );

  inline auto operator->() -> Cell *;
  inline auto operator*() -> Cell &;
  inline auto operator[]( size_t i ) -> Cell &;
};

//-----------------------------------------------------------------------------

inline CellIterator::CellIterator( Mesh & mesh )
  : MeshEntityIterator( mesh, mesh.topology_dimension() )
{
}

//-----------------------------------------------------------------------------

inline CellIterator::CellIterator( MeshEntity & entity )
  : MeshEntityIterator( entity, entity.mesh().topology_dimension() )
{
}

//-----------------------------------------------------------------------------

inline auto CellIterator::operator->() -> Cell *
{
  return static_cast< Cell * >( MeshEntityIterator::operator->() );
}

//-----------------------------------------------------------------------------

inline auto CellIterator::operator*() -> Cell &
{
  return *operator->();
}

//-----------------------------------------------------------------------------

inline auto CellIterator::operator[]( size_t i ) -> Cell &
{
  return static_cast< Cell & >( MeshEntityIterator::operator[]( i ) );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_CELL_ITERATOR_H */
