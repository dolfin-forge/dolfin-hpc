// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// Imported from licorne.
//
// First added:  2017-10-10
// Last changed: 2017-10-10

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

template<class E>
inline uint entity_dimension(Mesh& mesh);

//-----------------------------------------------------------------------------

template<class T, class E, uint N>
struct MeshValues : public MeshFunction<T>
{
  ///
  MeshValues(Mesh& mesh) :
    MeshFunction<T>(mesh, entity_dimension<E>(mesh))
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
  bool operator==(MeshValues<T, E, N> const& other)
  {
    return MeshFunction<T>::operator ==(other);
  }

  /// Equality
  bool operator!=(MeshValues<T, E, N> const& other)
  {
    return MeshFunction<T>::operator !=(other);
  }

  /// Return value size
  inline uint value_size() { return N; }

  /// Assignment operator
  MeshValues<T, E, N>& operator=(MeshValues<T, E, N> const& other)
  {
    MeshFunction<T>::operator=(other);
    return *this;
  }

  /// Assignment conversion operator
  template <class V>
  MeshValues<T, E, N>& operator=(MeshValues<V, E, N> const& other)
  {
    MeshFunction<T>::operator=(other);
    return *this;
  }

  /// Set all values to given value
  inline MeshValues<T, E, N>& operator=(T const& value)
  {
    MeshFunction<T>::operator=(value);
    return *this;
  }

  ///--- Value accessors

  // NOTE: operators below are defined for any value size but their use is
  //       not valid until MeshFunction supports vector values.

  /// Return value at given entity
  inline T& operator()(E& entity, uint i = 0)
  {
    dolfin_assert(this->values_);
    dolfin_assert(&entity.mesh() == this->mesh_);
    return this->values_[entity.index() * N + i];
  }

  /// Return value at given entity
  inline T const& operator()(E& entity, uint i = 0) const
  {
    dolfin_assert(this->values_);
    dolfin_assert(&entity.mesh() == this->mesh_);
    return this->values_[entity.index() * N + i];
  }

  /// Return value at given index
  inline T& operator()(uint index, uint i = 0)
  {
    dolfin_assert(this->values_);
    return this->values_[index * N + i];
  }

  /// Return value at given index
  inline T const& operator()(uint index, uint i = 0) const
  {
    dolfin_assert(this->values_);
    return this->values_[index * N + i];
  }

  ///--- Array accessors

  /// Return point to value at given entity
  inline T * operator[](E& entity)
  {
    dolfin_assert(this->values_);
    dolfin_assert(&entity.mesh() == this->mesh_);
    return this->values_ + entity.index() * N;
  }

  /// Return pointer to value at given entity
  inline T const * operator[](E& entity) const
  {
    dolfin_assert(this->values_);
    dolfin_assert(&entity.mesh() == this->mesh_);
    return this->values_ + entity.index() * N;
  }

  /// Return point to value at given index
  inline T * operator[](uint index)
  {
    dolfin_assert(this->values_);
    return this->values_ + index * N;
  }

  /// Return pointer to value at given index
  inline T const * operator[](uint index) const
  {
    dolfin_assert(this->values_);
    return this->values_ + index * N;
  }

};

//--- TEMPLATE SPECIALIZATIONS ------------------------------------------------

template<>
inline uint entity_dimension<Vertex>(Mesh& mesh)
{
  return 0;
}
template<>
inline uint entity_dimension<Edge>(Mesh& mesh)
{
  return 1;
}
template<>
inline uint entity_dimension<Face>(Mesh& mesh)
{
  return 2;
}
template<>
inline uint entity_dimension<Facet>(Mesh& mesh)
{
  return mesh.type().facet_dim();
}
template<>
inline uint entity_dimension<Cell>(Mesh& mesh)
{
  return mesh.type().dim();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_VALUES_H */

