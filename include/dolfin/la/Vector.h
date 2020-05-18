// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_VECTOR_H
#define __DOLFIN_VECTOR_H

#include <dolfin/la/GenericVector.h>

#include <dolfin/la/DefaultFactory.h>

namespace dolfin
{

/// This class provides the default DOLFIN vector class,
/// based on the default DOLFIN linear algebra backend.

class Vector : public GenericVector
{

public:

  /// Create empty vector
  Vector() :
      vector_(DefaultFactory::factory().createVector())
  {
  }

  /// Create vector of size N distributed by default
  explicit Vector(uint N) :
      vector_(DefaultFactory::factory().createVector())
  {
    vector_->init(N);
  }

  /// Create vector of size N distributed if specified
  explicit Vector(uint N, bool distributed) :
      vector_(DefaultFactory::factory().createVector())
  {
    vector_->init(N, distributed);
  }

  /// Copy constructor
  explicit Vector(Vector const& x) :
      vector_(x.vector_->copy())
  {
  }

  /// Destructor
  ~Vector()
  {
    delete vector_;
  }

  //--- Implementation of the GenericTensor interface ---

  /// Return copy of tensor
  Vector* copy() const
  {
    return new Vector(*this);
  }

  /// Set all entries to zero and keep any sparse structure
  void zero()
  {
    vector_->zero();
  }

  /// Finalize assembly of tensor
  void apply(FinalizeType finaltype = FINALIZE)
  {
    vector_->apply(finaltype);
  }

  /// Display tensor
  void disp(uint precision = 2) const
  {
    vector_->disp(precision);
  }

  //--- Implementation of the GenericVector interface ---

  /// Initialize vector of size N
  void init(uint N)
  {
    vector_->init(N);
  }

  /// Initialize vector of size N and distribute if specified
  void init(uint N, bool distributed)
  {
    vector_->init(N, distributed);
  }

  void init_ghosted(uint n, _ordered_set<uint>& indices,
                    _ordered_map<uint, uint>& map)
  {
    vector_->init_ghosted(n, indices, map);
  }

  /// Return size of vector
  uint size() const
  {
    return vector_->size();
  }

  /// Return local size of vector
  uint local_size() const
  {
    return vector_->local_size();
  }

  /// Return rank's offset into vector
  uint offset() const
  {
    return vector_->offset();
  }

  /// Get block of values
  void get(real* block, uint m, const uint* rows) const
  {
    vector_->get(block, m, rows);
  }

  /// Set block of values
  void set(const real* block, uint m, const uint* rows)
  {
    vector_->set(block, m, rows);
  }

  /// Add block of values
  void add(const real* block, uint m, const uint* rows)
  {
    vector_->add(block, m, rows);
  }

  /// Get all local values (not ghost entries)
  void get(real* values) const
  {
    vector_->get(values);
  }

  /// Set all local values (not ghost entries)
  void set(real* values)
  {
    vector_->set(values);
  }

  /// Add values to each local entry (not ghost entries)
  void add(real* values)
  {
    vector_->add(values);
  }

  /// Add multiple of given vector (AXPY operation)
  void axpy(real a, const GenericVector& x)
  {
    vector_->axpy(a, x);
  }

  /// Return inner product with given vector
  real inner(const GenericVector& x) const
  {
    return vector_->inner(x);
  }

  /// Return norm of vector
  real norm(VectorNormType type = l2) const
  {
    return vector_->norm(type);
  }

  /// Return minimum value of vector
  real min() const
  {
    return vector_->min();
  }

  /// Return maximum value of vector
  real max() const
  {
    return vector_->max();
  }

  /// Return pointwise operator op of vector and given vector x
  void pointwise(const GenericVector& x, VectorPointwiseOp op=pw_min) const
  {
    return vector_->pointwise(x, op);
  }

  /// Multiply vector by given number
  Vector& operator*=(real a)
  {
    *vector_ *= a;
    return *this;
  }

  /// Divide vector by given number
  Vector& operator/=(real a)
  {
    *this *= 1.0 / a;
    return *this;
  }

  /// Multiply vector by given vector component-wise
  Vector& operator*=(const GenericVector& x)
  {
    *vector_ *= x;
    return *this;
  }

  /// Add given vector
  Vector& operator+=(const GenericVector& x)
  {
    axpy(1.0, x);
    return *this;
  }

  /// Subtract given vector
  Vector& operator-=(const GenericVector& x)
  {
    axpy(-1.0, x);
    return *this;
  }

  /// Assignment operator
  GenericVector& operator=(const GenericVector& x)
  {
    *vector_ = x;
    return *this;
  }

  /// Assignment operator
  Vector& operator=(real a)
  {
    *vector_ = a;
    return *this;
  }

  //--- Special functions ---

  /// Return linear algebra backend factory
  LinearAlgebraFactory& factory() const
  {
    return vector_->factory();
  }

  //--- Special functions, intended for library use only ---

  /// Return concrete instance / unwrap (const)
  GenericVector const* instance() const
  {
    return vector_;
  }

  /// Return concrete instance / unwrap (non-const version)
  GenericVector* instance()
  {
    return vector_;
  }

  //--- Special Vector functions ---

  /// Assignment operator
  Vector& operator=(Vector const& x)
  {
    *vector_ = *x.vector_;
    return *this;
  }

private:

  // Pointer to concrete implementation
  GenericVector * const vector_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_VECTOR_H */
