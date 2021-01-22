// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_VALUES_H
#define __DOLFIN_MESH_VALUES_H

#include <dolfin/mesh/MeshFunction.h>

#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Edge.h>
#include <dolfin/mesh/entities/Face.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

/**
 *  @class  MeshValues
 *
 *  @brief  A MeshFunction defined given an entity (what it should have been).
 *
 */

template < class T, class E, size_t N >
struct MeshValues : public MeshFunction< T >
{
  ///
  MeshValues( Mesh & mesh, T val = static_cast< T >( 0 ) );

  ///
  template < class V >
  MeshValues( MeshValues< V, E, N > const & other );

  /// Equality
  auto operator==( MeshValues< T, E, N > const & other ) -> bool;

  /// Equality
  auto operator!=( MeshValues< T, E, N > const & other ) -> bool;

  /// Return value size
  auto value_size() -> size_t;

  /// Assignment operator
  auto operator=( MeshValues< T, E, N > const & other )
    -> MeshValues< T, E, N > &;

  /// Assignment conversion operator
  template < class V >
  auto operator=( MeshValues< V, E, N > const & other )
    -> MeshValues< T, E, N > &;

  /// Set all values to given value
  inline auto operator=( T const & value ) -> MeshValues< T, E, N > &;

  /// Swap operator
  template < class T_, class E_, size_t N_ >
  friend auto swap( MeshValues< T_, E_, N_ > & a, MeshValues< T_, E_, N_ > & b )
    -> void;

  ///--- Value accessors

  // NOTE: operators below are defined for any value size but their use is
  //       not valid until MeshFunction supports vector values.

  /// Return value at given entity
  inline auto operator()( E const & entity, size_t i = 0 ) -> T &;

  /// Return value at given entity
  inline auto operator()( E const & entity, size_t i = 0 ) const -> T const &;

  /// Return value at given index;
  inline auto operator()( size_t index, size_t i = 0 ) -> T &;

  /// Return value at given index
  inline auto operator()( size_t index, size_t i = 0 ) const -> T const &;

  ///--- std::vector accessors

  /// Return point to value at given entity
  inline auto operator[]( E & entity ) -> T *;

  /// Return pointer to value at given entity
  inline auto operator[]( E & entity ) const -> T const *;

  /// Return point to value at given index
  inline auto operator[]( size_t index ) -> T *;

  /// Return pointer to value at given index
  inline auto operator[]( size_t index ) const -> T const *;
};

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
MeshValues< T, E, N >::MeshValues( Mesh & mesh, T val )
  : MeshFunction< T >( mesh, entity_dimension< E >( mesh ), val )
{
  if ( N > 1 )
  {
    error( "MeshValues : vector values are unsupported for now." );
  }
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
template < class V >
MeshValues< T, E, N >::MeshValues( MeshValues< V, E, N > const & other )
  : MeshFunction< T >( other )
{
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto
  MeshValues< T, E, N >::operator==( MeshValues< T, E, N > const & other )
    -> bool
{
  return MeshFunction< T >::operator==( other );
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto
  MeshValues< T, E, N >::operator!=( MeshValues< T, E, N > const & other )
    -> bool
{
  return MeshFunction< T >::operator!=( other );
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto MeshValues< T, E, N >::value_size() -> size_t
{
  return N;
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto
  MeshValues< T, E, N >::operator=( MeshValues< T, E, N > const & other )
    -> MeshValues< T, E, N > &
{
  MeshFunction< T >::operator=( other );
  return *this;
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
template < class V >
inline auto
  MeshValues< T, E, N >::operator=( MeshValues< V, E, N > const & other )
    -> MeshValues< T, E, N > &
{
  MeshFunction< T >::operator=( other );
  return *this;
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto MeshValues< T, E, N >::operator=( T const & value )
  -> MeshValues< T, E, N > &
{
  MeshFunction< T >::operator=( value );
  return *this;
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto swap( MeshValues< T, E, N > & a, MeshValues< T, E, N > & b ) -> void
{
  swap( static_cast< MeshFunction< T > & >( a ),
        static_cast< MeshFunction< T > & >( b ) );
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto MeshValues< T, E, N >::operator()( E const & entity, size_t i )
  -> T &
{
  dolfin_assert( this->values_ );
  dolfin_assert( &entity.mesh() == this->mesh_ );
  dolfin_assert( entity.index() < this->size_ );
  return this->values_[entity.index() * N + i];
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto MeshValues< T, E, N >::operator()( E const & entity,
                                               size_t    i ) const -> T const &
{
  dolfin_assert( this->values_ );
  dolfin_assert( &entity.mesh() == this->mesh_ );
  dolfin_assert( entity.index() < this->size_ );
  return this->values_[entity.index() * N + i];
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto MeshValues< T, E, N >::operator()( size_t index, size_t i ) -> T &
{
  dolfin_assert( this->values_ );
  dolfin_assert( index < this->size_ );
  return this->values_[index * N + i];
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto MeshValues< T, E, N >::operator()( size_t index, size_t i ) const
  -> T const &
{
  dolfin_assert( this->values_ );
  dolfin_assert( index < this->size_ );
  return this->values_[index * N + i];
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto MeshValues< T, E, N >::operator[]( E & entity ) -> T *
{
  dolfin_assert( this->values_ );
  dolfin_assert( &entity.mesh() == this->mesh_ );
  dolfin_assert( entity.index() < this->size_ );
  return this->values_ + entity.index() * N;
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto MeshValues< T, E, N >::operator[]( E & entity ) const -> T const *
{
  dolfin_assert( this->values_ );
  dolfin_assert( &entity.mesh() == this->mesh_ );
  dolfin_assert( entity.index() < this->size_ );
  return this->values_ + entity.index() * N;
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto MeshValues< T, E, N >::operator[]( size_t index ) -> T *
{
  dolfin_assert( this->values_ );
  dolfin_assert( index < this->size_ );
  return this->values_ + index * N;
}

//-----------------------------------------------------------------------------

template < class T, class E, size_t N >
inline auto MeshValues< T, E, N >::operator[]( size_t index ) const -> T const *
{
  dolfin_assert( this->values_ );
  dolfin_assert( index < this->size_ );
  return this->values_ + index * N;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_VALUES_H */
