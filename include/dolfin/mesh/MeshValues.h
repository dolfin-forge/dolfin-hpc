// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_VALUES_H
#define __DOLFIN_MESH_VALUES_H

#include <dolfin/mesh/MeshFunction.h>

#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Face.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Cell.h>

namespace dolfin
{

/**
 *  @class  MeshValues
 *
 *  @brief  A MeshFunction defined given an entity (what it should have been).
 *
 */

//-----------------------------------------------------------------------------

template<class T, class E, uint N>
struct MeshValues : public MeshFunction<T>
{
  ///
  MeshValues(Mesh& mesh, T val = static_cast<T>(0)) :
    MeshFunction<T>(mesh, entity_dimension<E>(mesh), val)
  {
    if (N > 1)
    {
      error("MeshValues : vector values are unsupported for now.");
    }
  }

  ///
  template<class V>
  MeshValues(MeshValues<V, E, N> const& other) :
      MeshFunction<T>(other)
  {
  }

  /// Equality
  auto operator==(MeshValues<T, E, N> const& other) -> bool;

  /// Equality
  auto operator!=(MeshValues<T, E, N> const& other) -> bool;

  /// Return value size
  inline auto value_size() -> uint;

  /// Assignment operator
  auto operator=(MeshValues<T, E, N> const& other) -> MeshValues<T, E, N>&;

  /// Assignment conversion operator
  template <class V>
  auto operator=(MeshValues<V, E, N> const& other) -> MeshValues<T, E, N>&;

  /// Set all values to given value
  inline auto operator=(T const& value) -> MeshValues<T, E, N>&;

  /// Swap operator
  friend void swap( MeshValues<T,E,N>& a, MeshValues<T,E,N>& b )
  {
    swap( static_cast<MeshFunction<T>&>(a), static_cast<MeshFunction<T>&>(b) );
  }

  ///--- Value accessors

  // NOTE: operators below are defined for any value size but their use is
  //       not valid until MeshFunction supports vector values.

  /// Return value at given entity
  inline auto operator()(E const& entity, uint i = 0) -> T&;

  /// Return value at given entity
  inline auto operator()(E const& entity, uint i = 0) const -> T const&;

  /// Return value at given index;
  inline auto operator()(uint index, uint i = 0) -> T&;

  /// Return value at given index
  inline auto operator()(uint index, uint i = 0) const -> T const&;

  ///--- Array accessors

  /// Return point to value at given entity
  inline auto operator[](E& entity) -> T *;

  /// Return pointer to value at given entity
  inline auto operator[](E& entity) const -> T const *;

  /// Return point to value at given index
  inline auto operator[](uint index) -> T *;

  /// Return pointer to value at given index
  inline auto operator[](uint index) const -> T const *;

};

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto
  MeshValues< T, E, N >::operator==( MeshValues< T, E, N > const & other ) -> bool
{
  return MeshFunction< T >::operator==( other );
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto
  MeshValues< T, E, N >::operator!=( MeshValues< T, E, N > const & other ) -> bool
{
  return MeshFunction< T >::operator!=( other );
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto MeshValues< T, E, N >::value_size() -> uint
{
  return N;
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto
  MeshValues< T, E, N >::operator=( MeshValues< T, E, N > const & other ) -> MeshValues< T, E, N > &
{
  MeshFunction< T >::operator=( other );
  return *this;
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
template < class V >
inline auto
  MeshValues< T, E, N >::operator=( MeshValues< V, E, N > const & other ) -> MeshValues< T, E, N > &
{
  MeshFunction< T >::operator=( other );
  return *this;
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto
  MeshValues< T, E, N >::operator=( T const & value ) -> MeshValues< T, E, N > &
{
  MeshFunction< T >::operator=( value );
  return *this;
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto MeshValues< T, E, N >::operator()( E const & entity, uint i ) -> T &
{
  dolfin_assert( this->values_ );
  dolfin_assert( &entity.mesh() == this->mesh_ );
  dolfin_assert( entity.index() < this->size_ );
  return this->values_[entity.index() * N + i];
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto MeshValues< T, E, N >::operator()( E const & entity,
                                                    uint      i ) const -> T const &
{
  dolfin_assert( this->values_ );
  dolfin_assert( &entity.mesh() == this->mesh_ );
  dolfin_assert( entity.index() < this->size_ );
  return this->values_[entity.index() * N + i];
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto MeshValues< T, E, N >::operator()( uint index, uint i ) -> T &
{
  dolfin_assert( this->values_ );
  dolfin_assert( index < this->size_ );
  return this->values_[index * N + i];
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto MeshValues< T, E, N >::operator()( uint index, uint i ) const -> T const &
{
  dolfin_assert( this->values_ );
  dolfin_assert( index < this->size_ );
  return this->values_[index * N + i];
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto MeshValues< T, E, N >::operator[]( E & entity ) -> T *
{
  dolfin_assert( this->values_ );
  dolfin_assert( &entity.mesh() == this->mesh_ );
  dolfin_assert( entity.index() < this->size_ );
  return this->values_ + entity.index() * N;
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto MeshValues< T, E, N >::operator[]( E & entity ) const -> T const *
{
  dolfin_assert( this->values_ );
  dolfin_assert( &entity.mesh() == this->mesh_ );
  dolfin_assert( entity.index() < this->size_ );
  return this->values_ + entity.index() * N;
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto MeshValues< T, E, N >::operator[]( uint index ) -> T *
{
  dolfin_assert( this->values_ );
  dolfin_assert( index < this->size_ );
  return this->values_ + index * N;
}

//-----------------------------------------------------------------------------
template < class T, class E, uint N >
inline auto MeshValues< T, E, N >::operator[]( uint index ) const -> T const *
{
  dolfin_assert( this->values_ );
  dolfin_assert( index < this->size_ );
  return this->values_ + index * N;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_VALUES_H */

