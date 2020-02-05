// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

// This class suffered from some cleanup and was cured from segmentation faults
// on initialization to zero size, which in the course of 2006-2012 was the
// cause of many ugly workarounds as well as erratic behaviour of parallel
// algorithms.

#ifndef __DOLFIN_MESH_FUNCTION_H
#define __DOLFIN_MESH_FUNCTION_H

#include <dolfin/common/types.h>
#include <dolfin/io/File.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEntity.h>

#include <algorithm>

namespace dolfin
{

/**
 *  @class  MeshFunction
 *
 *  @brief  A MeshFunction is a function that can be evaluated at a set of
 *          mesh entities.
 *          A MeshFunction is discrete and is only defined at the set of mesh
 *          entities of a fixed topological dimension.
 *          A MeshFunction may for example be used to store a global numbering
 *          scheme for the entities of a (parallel) mesh, marking sub domains or
 *          boolean markers for mesh refinement.
 *
 *          MeshFunctions are default initialized to the zero value of the
 *          underlying type:
 *            - bool : false
 *            - uint : 0
 *            - real : 0.0
 */

template<typename T>
class XMLMeshFunction;

template<class T>
class MeshFunction
{

public:

  /// Copy constructor
  explicit
  MeshFunction(MeshFunction<T> const& other) :
      mesh_(0),
      dim_(0),
      size_(0),
      values_(NULL)
  {
    MeshFunction<T>::operator=(other);
  }

  /// Copy constructor
  template <class V>
  MeshFunction(MeshFunction<V> const& other) :
      mesh_(0),
      dim_(0),
      size_(0),
      values_(NULL)
  {
    MeshFunction<T>::operator=(other);
  }

  /// Destructor
  virtual ~MeshFunction()
  {
    delete[] values_;
  }

  ///
  bool empty() const
  {
    return (values_ == NULL);
  }

  /// Assignment operator
  MeshFunction<T>& operator=(MeshFunction<T> const& other)
  {
    if(this != &other)
    {
      if(this->mesh_ != other.mesh_ || this->size_ != other.size_ )
      {
        init(other.mesh_, other.dim_, other.size_);
      }
      this->dim_ = other.dim_;
      std::copy(other.values_, other.values_ + size_, values_);
    }
    return *this;
  }

  /// Required for assignment operator
  friend class MeshFunction<bool>;
  friend class XMLMeshFunction<bool>;

  friend class MeshFunction<int>;
  friend class XMLMeshFunction<int>;

  friend class MeshFunction<uint>;
  friend class XMLMeshFunction<uint>;

  friend class MeshFunction<float>;
  friend class XMLMeshFunction<float>;

  friend class MeshFunction<real>;
  friend class XMLMeshFunction<real>;

  /// Assignment conversion operator
  template <class V>
  MeshFunction<T>& operator=(MeshFunction<V> const& other)
  {
    if(this->mesh_ != other.mesh_ || this->size_ != other.size_ )
    {
      init(other.mesh_, other.dim_, other.size_);
    }
    this->dim_ = other.dim_;
    std::transform(other.values_, other.values_ + size_, values_, cast<V>() );
    return *this;
  }

  /// Equality
  bool operator==(MeshFunction<T> const& other)
  {
    if(this == &other)
    {
      return true;
    }
    if (dim_ != other.dim_)
    {
      return false;
    }
    if (size_ != other.size_)
    {
      return false;
    }
    if (size_ == 0)
    {
      return true;
    }
    for (uint i = 0; i < size_; ++i)
    {
      if (values_[i] != other.values_[i])
      {
        return false;
      }
    }
    return true;
  }

  /// Equality with cast
  template <class V>
  bool operator==(MeshFunction<V> const& other)
  {
    if(this == &other)
    {
      return true;
    }
    if (dim_ != other.dim_)
    {
      return false;
    }
    if (size_ != other.size_)
    {
      return false;
    }
    if (size_ == 0)
    {
      return true;
    }
    for (uint i = 0; i < size_; ++i)
    {
      if (values_[i] != cast<V>(other.values_[i]))
      {
        return false;
      }
    }
    return true;
  }

  /// Equality
  bool operator!=(MeshFunction<T> const& other)
  {
    return !(*this == other);
  }

  /// Return mesh associated with mesh function
  inline Mesh& mesh() const
  {
    dolfin_assert(mesh_);
    return *mesh_;
  }

  /// Return topological dimension
  inline uint dim() const
  {
    return dim_;
  }

  /// Return size (number of entities)
  inline uint size() const
  {
    return size_;
  }

  /// Access value at given index
  inline T& operator()(uidx index)
  {
    dolfin_assert(values_);
    dolfin_assert(index < size_);
    return values_[index];
  }

  /// Access value at given index (const)
  inline T const& operator()(uidx index) const
  {
    dolfin_assert(values_);
    dolfin_assert(index < size_);
    return values_[index];
  }

  /// Return value at given entity
  inline T& operator()(MeshEntity& entity)
  {
    dolfin_assert(values_);
    dolfin_assert(&entity.mesh() == mesh_);
    dolfin_assert(entity.dim() == dim_);
    dolfin_assert(entity.index() < size_);
    return values_[entity.index()];
  }

  /// Return value at given entity
  inline T const& operator()(MeshEntity& entity) const
  {
    dolfin_assert(values_);
    dolfin_assert(&entity.mesh() == mesh_);
    dolfin_assert(entity.dim() == dim_);
    dolfin_assert(entity.index() < size_);
    return values_[entity.index()];
  }

  /// Set all values to given value
  MeshFunction<T>& operator=(T const& value)
  {
    dolfin_assert(!((values_ == NULL) && (size_>0)));
    std::fill_n(values_, size_, value);
    return *this;
  }

  /// Swap instances
  friend void swap( MeshFunction<T>& a, MeshFunction<T>& b )
  {
    using std::swap;

		swap( a.mesh_,   b.mesh_   );
		swap( a.dim_,    b.dim_    );
		swap( a.size_,   b.size_   );
		swap( a.values_, b.values_ );
	}

  /// Display mesh function data
  void disp() const
  {
    section("MeshFunction");
    cout << "Topological dimension: " << dim_ << "\n";
    cout << "Number of values     : " << size_ << "\n";
    cout << "\n";
    for (uint i = 0; i < size_; ++i)
    {
      cout << "(" << dim_ << ", " << i << "): " << values_[i] << "\n";
    }
    end();
  }

protected:

  /// Create scalar mesh function on given mesh of given dimension
  MeshFunction(Mesh& mesh, uint dim, T val = static_cast<T>(0)) :
      mesh_(0),
      dim_(0),
      size_(0),
      values_(NULL)
  {
    init(&mesh, dim, mesh.size(dim));
    std::fill_n(values_, size_, val);
  }

  /// Initialize mesh function for given topological dimension of given size
  void init(Mesh * mesh, uint dim, uint size)
  {
    mesh_ = mesh;
    dim_ = dim;
    size_ = size;
    delete[] values_;
    values_ = NULL;
    if(size_ > 0)
    {
      values_ = new T[size];
    }
  }

  /// Cast operators
  template<class V>
  struct cast
  {
    inline T operator()(V x) const { return static_cast<T>(x); }
  };

  /// The mesh
  Mesh * mesh_;

  /// Topological dimension
  uint dim_;

  /// Number of mesh entities
  uint size_;

  /// Values at the set of mesh entities
  T * values_;

};

//--- TEMPLATE SPECIALIZATIONS ------------------------------------------------

// Copyright (C) 2013 Balthasar Reuter.
// Licensed under the GNU LGPL Version 2.1.
/// Helper function that performs symmetric rounding to closest integer

template<> template<> inline
bool MeshFunction<bool>::cast<float>::operator()(float x) const
{
  return static_cast<bool>(x > 0);
}

template<> template<> inline
bool MeshFunction<bool>::cast<real>::operator()(real x) const
{
  return static_cast<bool>(x > 0);
}

template<> template<> inline
int MeshFunction<int>::cast<float>::operator()(float x) const
{
  return static_cast<int>((x > 0) ? std::floor(x + 0.5) : std::ceil(x - 0.5));
}

template<> template<> inline
int MeshFunction<int>::cast<real>::operator()(real x) const
{
  return static_cast<int>((x > 0) ? std::floor(x + 0.5) : std::ceil(x - 0.5));
}

template<> template<> inline
uint MeshFunction<uint>::cast<float>::operator()(float x) const
{
  return static_cast<uint>(std::floor(x + 0.5));
}

template<> template<> inline
uint MeshFunction<uint>::cast<real>::operator()(real x) const
{
  return static_cast<uint>(std::floor(x + 0.5));
}

template<> template<> inline
uint MeshFunction<uint>::cast<int>::operator()(int x) const
{
  return static_cast<uint>((x > 0) ? x : 0);
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_FUNCTION_H */
