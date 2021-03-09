// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FACE_IITERATOR_H
#define __DOLFIN_FACE_IITERATOR_H

#include <dolfin/mesh/entities/Face.h>
#include <dolfin/mesh/entities/MeshEntity.h>
#include <dolfin/mesh/entities/iterators/MeshEntityIterator.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

/// A FaceIterator is a MeshEntityIterator of topological dimension 2.

class FaceIterator : public MeshEntityIterator
{
public:
  FaceIterator( Mesh & mesh );
  FaceIterator( MeshEntity & entity );

  inline auto operator*() -> Face &;
  inline auto operator->() -> Face *;
};

//-----------------------------------------------------------------------------

inline FaceIterator::FaceIterator( Mesh & mesh )
  : MeshEntityIterator( mesh, 2 )
{
}

//-----------------------------------------------------------------------------

inline FaceIterator::FaceIterator( MeshEntity & entity )
  : MeshEntityIterator( entity, 2 )
{
}

//-----------------------------------------------------------------------------

inline auto FaceIterator::operator*() -> Face &
{
  return *operator->();
}

//-----------------------------------------------------------------------------

inline auto FaceIterator::operator->() -> Face *
{
  return static_cast< Face * >( MeshEntityIterator::operator->() );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_FACE_IITERATOR_H */
