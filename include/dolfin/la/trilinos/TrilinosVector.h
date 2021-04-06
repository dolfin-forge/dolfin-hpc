// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TRILINOS_VECTOR_H
#define __DOLFIN_TRILINOS_VECTOR_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/common/Variable.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/trilinos/TrilinosObject.h>

#include <Tpetra_MultiVector.hpp>

namespace dolfin
{

namespace trilinos
{

class Matrix;

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
  /// local index type
  using LO = int;

  /// global index type
  using GO = size_t;

  /// array indexing starts at 0
  static constexpr GO const indexBase = 0;

  /// Tpetra node type
  using TPNode      = Tpetra::MultiVector<>::node_type;

  /// TpetraVector map type (local index, global index)
  using TPMap       = Tpetra::Map< LO, GO, TPNode >;

  /// TpetraVector vector type (scalar, local index, global index, node)
  using TPVector    = Tpetra::MultiVector< real, LO, GO, TPNode >;

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
  auto copy() const -> Vector * override;

  /// Set all entries to zero and keep any sparse structure
  auto zero() -> void override;

  /// Finalize assembly of tensor
  auto apply( FinalizeType ) -> void override;

  /// Display tensor
  auto disp( size_t precision = 0 ) const -> void override;

  //--- Implementation of the GenericVector interface ---

  /// Return size of vector
  auto size() const -> size_t override;

  /// Return local size of vector
  auto local_size() const -> size_t override;

  /// Return rank's offset into vector
  auto offset() const -> size_t override;

  /// Initialize vector of local size N, distributed by default
  auto init( size_t N ) -> void override;

  /// Initialize vector of local size N, distributed if specified
  auto init( size_t N, bool distributed ) -> void override;

  /// Initialize ghost entries
  auto init_ghosted( size_t                           n,
                     _ordered_set< size_t > &         indices,
                     _ordered_map< size_t, size_t > & map ) -> void override;

  /// Get block of values
  auto get( real * block, size_t m, const size_t * rows ) const -> void override;

  /// Set block of values
  auto set( const real * block, size_t m, const size_t * rows ) -> void override;

  /// Add block of values
  auto add( const real * block, size_t m, const size_t * rows ) -> void override;

  /// Get all values
  auto get( real * values ) const -> void override;

  /// Set all values
  auto set( real * values ) -> void override;

  /// Add values to each entry
  auto add( real * values ) -> void override;

  /// Add multiple of given vector (AXPY operation)
  auto axpy( real a, const GenericVector & x ) -> void override;

  /// Add multiple of given vector (y=a*x+b*y)
  auto axpby( real a, const GenericVector & x, real b ) -> void override;

  /// Add multiple of given vector (w=a*x+y)
  auto waxpy( real a, const GenericVector & x, const GenericVector & y ) -> void override;

  /// Add multiple of given vector (z=a*x+b*y+c*z)
  auto axpbypcz( real a, const GenericVector & x,
                 real b, const GenericVector & y,
                 real c ) -> void override;

  /// Return inner product with given vector
  auto inner( const GenericVector & v ) const -> real override;

  /// Return norm of vector
  auto norm( VectorNormType type = l2 ) const -> real override;

  /// Return minimum value of vector
  auto min() const -> real override;

  /// Return maximum value of vector
  auto max() const -> real override;

  /// Return pointwise operator op of vector and given vector x
  auto pointwise( const GenericVector & x,
                  VectorPointwiseOp     op = pw_min ) const -> void override;

  /// Multiply vector by given number
  auto operator*=( real a ) -> Vector & override;

  /// Divide vector by given number
  auto operator/=( real a ) -> Vector & override;

  /// Multiply vector by given vector component-wise
  auto operator*=( const GenericVector & x ) -> Vector & override;

  /// Add given vector
  auto operator+=( const GenericVector & x ) -> Vector & override;

  /// Subtract given vector
  auto operator-=( const GenericVector & x ) -> Vector & override;

  /// Assignment operator
  auto operator=( Vector const & x ) -> Vector &;

  /// Assignment operator
  auto operator=( GenericVector const & x ) -> Vector & override;

  /// Assignment operator
  auto operator=( real a ) -> Vector & override;

  /// Return Trilinos Vec pointer
  auto vec() const -> TPVectorPtr;

  //--- Special functions ---

  /// Return linear algebra backend factory
  auto factory() const -> LinearAlgebraFactory & override;

private:
  friend trilinos::Matrix;

  //
  auto clear() -> void;

  // Tpetra multivector
  TPVectorPtr vec_;
};

} // end namespace trilinos

} // end namespace dolfin

#endif // HAVE_TRILINOS

#endif // __DOLFIN_TRILINOS_VECTOR_H
