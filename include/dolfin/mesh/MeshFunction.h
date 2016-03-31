// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman, 2007.
//
// First added:  2006-05-22
// Last changed: 2008-05-21

#ifndef __DOLFIN_MESH_FUNCTION_H
#define __DOLFIN_MESH_FUNCTION_H

#include <dolfin/common/types.h>
#include <dolfin/io/File.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEntity.h>

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

template<class T>
class MeshFunction
{
public:

  /// Create empty mesh function
  MeshFunction() :
      values_(NULL),
      mesh_(0),
      dim_(0),
      size_(0)
  {
  }

  /// Create empty mesh function on given mesh
  MeshFunction(Mesh& mesh) :
      values_(NULL),
      mesh_(&mesh),
      dim_(0),
      size_(0)
  {
  }

  /// Create mesh function on given mesh of given dimension
  MeshFunction(Mesh& mesh, uint dim) :
      values_(NULL),
      mesh_(&mesh),
      dim_(0),
      size_(0)
  {
    init(dim);
  }

  /// Create function from data file
  MeshFunction(Mesh& mesh, const std::string filename) :
      values_(NULL),
      mesh_(&mesh),
      dim_(0),
      size_(0)
  {
    File file(filename);
    file >> *this;
  }

  /// Destructor
  ~MeshFunction()
  {
    delete[] values_;
  }

  /// Return mesh associated with mesh function
  inline Mesh& mesh()
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

  /// Return array of values
  inline const T* values() const
  {
    return values_;
  }

  /// Return array of values
  inline T* values()
  {
    return values_;
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
  inline const T& operator()(MeshEntity& entity) const
  {
    dolfin_assert(values_);
    dolfin_assert(&entity.mesh() == mesh_);
    dolfin_assert(entity.dim() == dim_);
    dolfin_assert(entity.index() < size_);
    return values_[entity.index()];
  }

  /// Set all values to given value
  const MeshFunction<T>& operator=(const T& value)
  {
    dolfin_assert(values_);
    for (uint i = 0; i < size_; ++i)
    {
      values_[i] = value;
    }
    return *this;
  }

  /// Equality
  bool operator==(const MeshFunction<T>& other)
  {
    if(this == &other)
    {
      return true;
    }
    if (size_ != other.size_)
    {
      return false;
    }
    if (size_ == 0)
    {
      return true;
    }
    bool cmp = true;
    for (uint i = 0; i < size_; i++)
    {
      cmp &= (values_[i] == other.values_[i]);
    }
    return cmp;
  }

  /// Equality
  bool operator!=(const MeshFunction<T>& other)
  {
    return !(*this == other);
  }

  /// Initialize mesh function for given topological dimension
  void init(uint dim)
  {
    if (!mesh_)
    {
      error("Mesh undefined, unable to initialize mesh function.");
    }
    mesh_->init(dim);
    init(*mesh_, dim, mesh_->size(dim));
  }

  /// Initialize mesh function for given topological dimension of given size
  void init(uint dim, uint size)
  {
    if (!mesh_)
    {
      error("Mesh undefined, unable to initialize mesh function.");
    }
    mesh_->init(dim);
    init(*mesh_, dim, size);
  }

  /// Initialize mesh function for given topological dimension
  void init(Mesh& mesh, uint dim)
  {
    mesh.init(dim);
    init(mesh, dim, mesh.size(dim));
  }

  /// Initialize mesh function for given topological dimension of given size
  void init(Mesh& mesh, uint dim, uint size)
  {
    // Initialize mesh for entities of given dimension
    mesh.init(dim);
    dolfin_assert(mesh.size(dim) == size);

    // Initialize data
    mesh_ = &mesh;
    dim_ = dim;
    size_ = size;
    delete[] values_;
    values_ = NULL;

    if(size_ > 0)
    {
      values_ = new T[size];
      std::fill(values_, values_ + size, static_cast<T>(0));
    }
  }

  /// Get value at given entity
  inline T get(const MeshEntity& entity) const
  {
    dolfin_assert(values_);
    dolfin_assert(&entity.mesh() == mesh_);
    dolfin_assert(entity.dim() == dim_);
    dolfin_assert(entity.index() < size_);
    return values_[entity.index()];
  }

  /// Get value at given entity
  inline T get(uint index) const
  {
    dolfin_assert(values_);
    dolfin_assert(index < size_);
    return values_[index];
  }

  /// Set value at given entity
  inline void set(const MeshEntity& entity, const T& value)
  {
    dolfin_assert(values_);
    dolfin_assert(&entity.mesh() == mesh_);
    dolfin_assert(entity.dim() == dim_);
    dolfin_assert(entity.index() < size_);
    values_[entity.index()] = value;
  }

  /// Set value at given entity
  inline void set(uint index, const T& value)
  {
    dolfin_assert(values_);
    dolfin_assert(index < size_);
    values_[index] = value;
  }

  /// Display mesh function data
  void disp() const
  {
    section("Mesh function data");
    cout << "Topological dimension: " << dim_ << endl;
    cout << "Number of values:      " << size_ << endl;
    cout << endl;
    for (uint i = 0; i < size_; ++i)
    {
      cout << "(" << dim_ << ", " << i << "): " << values_[i] << endl;
    }
    end();
  }

private:

  /// Values at the set of mesh entities
  T * values_;

  /// The mesh
  Mesh* mesh_;

  /// Topological dimension
  uint dim_;

  /// Number of mesh entities
  uint size_;
};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_FUNCTION_H */
