// Copyright (C) 2020 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TRILINOS_VECTOR_H
#define __DOLFIN_TRILINOS_VECTOR_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/common/Variable.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/trilinos/TrilinosObject.h>

#include <Tpetra_MultiVector_def.hpp>

namespace dolfin
{

namespace trilinos
{

/// This class provides a simple vector class based on Trilinos/Tpetra.
/// It is a simple wrapper for a Tpetra vector pointer (Vec)
/// implementing the GenericVector interface.
///
/// The interface is intentionally simple. For advanced usage,
/// access the Trilinos/Tpetra Vec pointer using the function vec() and
/// use the standard Trilinos/Tpetra interface.

class Vector : public GenericVector, public Object, public Variable
{
public:
  /// TpetraVector map type (local index, global index)
  using TPMap       = Tpetra::Map< int, int, Tpetra::MultiVector<>::node_type >;

  /// TpetraVector vector type (scalar, local index, global index, node)
  using TPVector    = Tpetra::MultiVector< real, int, int, Tpetra::MultiVector<>::node_type >;

  /// smart pointer to a TPVector
  using TPVectorPtr = Teuchos::RCP< TPVector >;

public:
  /// Create empty vector
  Vector();

  /// Create vector of local size N
  explicit Vector( size_t N, bool distributed = true );

  /// Copy constructor
  explicit Vector( Vector const & copy );

  /// Create vector from given Tpetra Vec pointer
  // explicit Vector( TPVector x );

  /// Destructor
  ~Vector() override;

  //--- Implementation of the GenericTensor interface ---

  /// Return copy of tensor
  Vector * copy() const override;

  /// Set all entries to zero and keep any sparse structure
  void zero() override;

  /// Finalize assembly of tensor
  void apply( FinalizeType ) override;

  /// Display tensor
  void disp( size_t precision = 0 ) const override;

  //--- Implementation of the GenericVector interface ---

  /// Return size of vector
  size_t size() const override;

  /// Return local size of vector
  size_t local_size() const override;

  /// Return rank's offset into vector
  size_t offset() const override;

  /// Initialize vector of local size N, distributed by default
  void init( size_t N ) override;

  /// Initialize vector of local size N, distributed if specified
  void init( size_t N, bool distributed ) override;

  /// Initialize ghost entries
  void init_ghosted( size_t                           n,
                     _ordered_set< size_t > &         indices,
                     _ordered_map< size_t, size_t > & map ) override;

  /// Get block of values
  void get( real * block, size_t m, const size_t * rows ) const override;

  /// Set block of values
  void set( const real * block, size_t m, const size_t * rows ) override;

  /// Add block of values
  void add( const real * block, size_t m, const size_t * rows ) override;

  /// Get all values
  void get( real * values ) const override;

  /// Set all values
  void set( real * values ) override;

  /// Add values to each entry
  void add( real * values ) override;

  /// Add multiple of given vector (AXPY operation)
  void axpy( real a, const GenericVector & x ) override;

  /// Add multiple of given vector (y=a*x+b*y)
  void axpby( real a, const GenericVector & x, real b ) override;

  /// Add multiple of given vector (w=a*x+y)
  void waxpy( real a, const GenericVector & x, const GenericVector & y ) override;

  /// Add multiple of given vector (z=a*x+b*y+c*z)
  void axpbypcz( real a, const GenericVector & x,
                 real b, const GenericVector & y,
                 real c ) override;

  /// Return inner product with given vector
  real inner( const GenericVector & v ) const override;

  /// Return norm of vector
  real norm( VectorNormType type = l2 ) const override;

  /// Return minimum value of vector
  real min() const override;

  /// Return maximum value of vector
  real max() const override;

  /// Return pointwise operator op of vector and given vector x
  void pointwise( const GenericVector & x,
                  VectorPointwiseOp     op = pw_min ) const override;

  /// Multiply vector by given number
  Vector & operator*=( real a ) override;

  /// Divide vector by given number
  Vector & operator/=( real a ) override;

  /// Multiply vector by given vector component-wise
  Vector & operator*=( const GenericVector & x ) override;

  /// Add given vector
  Vector & operator+=( const GenericVector & x ) override;

  /// Subtract given vector
  Vector & operator-=( const GenericVector & x ) override;

  /// Assignment operator
  Vector & operator=( Vector const & x );

  /// Assignment operator
  Vector & operator=( GenericVector const & x ) override;

  /// Assignment operator
  Vector & operator=( real a ) override;

  /// Return Trilinos Vec pointer
  TPVectorPtr vec() const;

  //--- Special functions ---

  /// Return linear algebra backend factory
  LinearAlgebraFactory & factory() const override;

private:
  //
  void clear();

  // Tpetra multivector - actually a view into the ghosted vector, below
  TPVectorPtr x_;

  // Tpetra multivector with extra rows for ghost values
  TPVectorPtr x_ghosted_;
};

} // end namespace trilinos

} // end namespace dolfin

#endif // HAVE_TRILINOS

#endif // __DOLFIN_TRILINOS_VECTOR_H
