// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

// This class was initially cured in 2012 from a terrible disease "infinite
// loop on zero-size connectivity".
// The absence of robustness of the original class is fixed by "on-demand"
// computation of connectivities in the mesh topology.

#ifndef __DOLFIN_MESH_ENTITY_ITERATOR_H
#define __DOLFIN_MESH_ENTITY_ITERATOR_H

#include <dolfin/common/types.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/MeshEntity.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

/// MeshEntityIterator provides a common iterator for mesh entities
/// over meshes, boundaries and incidence relations. The basic use
/// is illustrated below.
///
/// The following example shows how to iterate over all mesh entities
/// of a mesh of topological dimension dim:
///
///     for (MeshEntityIterator e(mesh, dim); !e.end(); ++e)
///     {
///       e->foo();
///     }
///
/// The following example shows how to iterate over mesh entities of
/// topological dimension dim connected (incident) to some mesh entity f:
///
///     for (MeshEntityIterator e(f, dim); !e.end(); ++e)
///     {
///       e->foo();
///     }
///
/// In addition to the general iterator, a set of specific named iterators
/// are provided for entities of type Vertex, Edge, Face, Facet and Cell.
/// These iterators are defined along with their respective classes.

class MeshEntityIterator
{
public:
  /// Create iterator for mesh entities over given topological dimension
  MeshEntityIterator( Mesh & mesh, size_t dim );

  /// Create iterator for entities of given dimension connected to given entity
  MeshEntityIterator( MeshEntity & entity, size_t dim );

  /// Copy constructor is deleted to disallow usage. it would allow code like
  ///
  /// for (CellIterator c0(mesh); !c0.end(); ++c0)
  ///   for (CellIterator c1(c0); !c1.end(); ++c1)
  ///      ...
  ///
  /// c1 looks to be an iterator over the entities around c0 when it is in
  /// fact a copy of c0.
  MeshEntityIterator( MeshEntityIterator & other ) = delete;

  /// Destructor
  virtual ~MeshEntityIterator() = default;

  /// Step to next mesh entity (prefix increment)
  auto operator++() -> MeshEntityIterator &;

  /// Return current position
  auto pos() const -> size_t;

  /// Check if iterator has reached the end
  auto end() const -> bool;

  /// Member access operator
  auto operator->() -> MeshEntity *;

  /// Dereference operator
  auto operator*() -> MeshEntity &;

  /// Access operator at position relative to begin
  auto operator[]( size_t i ) -> MeshEntity &;

private:
  // Mesh entity
  MeshEntity entity_;

  // Current position
  size_t pos_;

  // End position
  size_t const end_;

  // Mapping from pos to index (if any)
  std::vector< size_t > const * index_;
};

//--- INLINES -----------------------------------------------------------------

inline MeshEntityIterator::MeshEntityIterator( Mesh & mesh, size_t dim )
  : entity_( mesh, dim, 0 )
  , pos_( 0 )
  , end_( mesh.size( dim ) )
  , index_( nullptr )
{
}

//-----------------------------------------------------------------------------
inline MeshEntityIterator::MeshEntityIterator( MeshEntity & entity, size_t dim )
  : entity_( entity.mesh(), dim, 0 )
  , pos_( 0 )
  , end_( entity.num_entities( dim ) )
  , index_( &entity.entities( dim ) )
{
}

//-----------------------------------------------------------------------------

inline auto MeshEntityIterator::operator++() -> MeshEntityIterator &
{
  ++pos_;
  return *this;
}

//-----------------------------------------------------------------------------

inline auto MeshEntityIterator::pos() const -> size_t
{
  return pos_;
}

//-----------------------------------------------------------------------------

inline auto MeshEntityIterator::end() const -> bool
{
  return pos_ >= end_;
}

//-----------------------------------------------------------------------------

inline auto MeshEntityIterator::operator->() -> MeshEntity *
{
  // WARNING: index is only updated if iterator is dereferenced
  entity_.index_ = ( index_ == nullptr ) ? pos_ : ( *index_ )[pos_];
  return &entity_;
}

//-----------------------------------------------------------------------------

inline auto MeshEntityIterator::operator*() -> MeshEntity &
{
  // WARNING: index is only updated if iterator is dereferenced
  entity_.index_ = ( index_ == nullptr ) ? pos_ : ( *index_ )[pos_];
  return entity_;
}

//-----------------------------------------------------------------------------

inline auto MeshEntityIterator::operator[]( size_t i ) -> MeshEntity &
{
  dolfin_assert( i < end_ );
  pos_           = i;
  entity_.index_ = ( index_ == nullptr ) ? pos_ : ( *index_ )[pos_];
  return entity_;
}

//-----------------------------------------------------------------------------

}

#endif
