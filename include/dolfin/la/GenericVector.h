// Copyright (C) 2006-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_GENERIC_VECTOR_H
#define __DOLFIN_GENERIC_VECTOR_H

#include <dolfin/la/GenericSparsityPattern.h>
#include <dolfin/la/GenericTensor.h>
#include <dolfin/la/VectorNormType.h>
#include <dolfin/la/VectorPointwiseOp.h>

#include <dolfin/common/maybe_unused.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

/// This class defines a common interface for vectors.

class GenericVector : public GenericTensor
{
public:
  /// Destructor
  ~GenericVector() override = default;

  //--- Implementation of the GenericTensor interface ---

  /// Initialize zero tensor using sparsity pattern
  inline void init( const GenericSparsityPattern & sparsity_pattern ) override;

  /// Return copy of tensor
  auto copy() const -> GenericVector * override = 0;

  /// Return tensor rank (number of dimensions)
  auto rank() const -> size_t override;

  /// Return size of given dimension
  inline auto size( size_t dim ) const -> size_t override;

  /// Get block of values
  inline void get( real *                 block,
                   const size_t *         num_rows,
                   const size_t * const * rows ) const override;

  /// Set block of values
  inline void set( const real *           block,
                   const size_t *         num_rows,
                   const size_t * const * rows ) override;

  /// Add block of values
  inline void add( const real *           block,
                   const size_t *         num_rows,
                   const size_t * const * rows ) override;

  /// Set all entries to zero and keep any sparse structure
  void zero() override = 0;

  /// Finalize assembly of tensor
  void apply( FinalizeType finaltype = FINALIZE ) override = 0;

  /// Display tensor
  void disp( size_t precision = 2 ) const override = 0;

  //--- Vector interface ---

  /// Initialize vector of size N
  virtual void init( size_t N ) = 0;

  /// Initialize vector of size N and distribute if specified
  virtual void init( size_t N, bool distributed ) = 0;

  /// Initialize ghost entries
  virtual void init_ghosted( size_t                           n,
                             _ordered_set< size_t > &         indices,
                             _ordered_map< size_t, size_t > & map ) = 0;

  /// Return size of vector
  virtual auto size() const -> size_t = 0;

  /// Return local size of vector
  virtual auto local_size() const -> size_t = 0;

  /// Return rank's offset into vector
  virtual auto offset() const -> size_t = 0;

  /// Get block of values
  virtual void get( real * block, size_t m, const size_t * rows ) const = 0;

  /// Set block of values
  virtual void set( const real * block, size_t m, const size_t * rows ) = 0;

  /// Add block of values
  virtual void add( const real * block, size_t m, const size_t * rows ) = 0;

  /// Get all values
  virtual void get( real * values ) const = 0;

  /// Set all values
  virtual void set( real * values ) = 0;

  /// Add values to each entry
  virtual void add( real * values ) = 0;

  /// Add multiple of given vector (y=a*x+y)
  virtual void axpy( real a, const GenericVector & x ) = 0;

  /// Add multiple of given vector (y=a*x+b*y)
  virtual void axpby( real a, const GenericVector & x, real b ) = 0;

  /// Add multiple of given vector (w=a*x+y)
  virtual void
    waxpy( real a, const GenericVector & x, const GenericVector & y ) = 0;

  /// Add multiple of given vector (z=a*x+b*y+c*z)
  virtual void axpbypcz( real                  a,
                         const GenericVector & x,
                         real                  b,
                         const GenericVector & y,
                         real                  c ) = 0;

  /// Return inner product with given vector
  virtual auto inner( const GenericVector & x ) const -> real = 0;

  /// Return norm of vector
  virtual auto norm( VectorNormType type = l2 ) const -> real = 0;

  /// Return minimum value of vector
  virtual auto min() const -> real = 0;

  /// Return maximum value of vector
  virtual auto max() const -> real = 0;

  /// Return pointwise operator op of vector and given vector x
  virtual void pointwise( const GenericVector & x,
                          VectorPointwiseOp     op = pw_min ) const = 0;

  /// Multiply vector by given number
  virtual auto operator*=( real a ) -> const GenericVector & = 0;

  /// Divide vector by given number
  virtual auto operator/=( real a ) -> const GenericVector & = 0;

  /// Multiply vector by given vector component-wise
  virtual auto operator*=( const GenericVector & x )
    -> const GenericVector & = 0;

  /// Add given vector
  virtual auto operator+=( const GenericVector & x )
    -> const GenericVector & = 0;

  /// Subtract given vector
  virtual auto operator-=( const GenericVector & x )
    -> const GenericVector & = 0;

  /// Assignment operator
  virtual auto operator      =( const GenericVector & x )
    -> const GenericVector & = 0;

  /// Assignment operator
  virtual auto operator=( real a ) -> const GenericVector & = 0;

  //--- Convenience functions ---

  /// Get value of given entry
  virtual auto operator[]( size_t i ) const -> real;

  /// Get value of given entry
  virtual auto getitem( size_t i ) const -> real;

  /// Set given entry to value
  virtual void setitem( size_t i, real value );
};

//-----------------------------------------------------------------------------

inline void GenericVector::init( const GenericSparsityPattern & sparsity_pattern )
{
  init( sparsity_pattern.size( 0 ), sparsity_pattern.is_distributed() );
}

//-----------------------------------------------------------------------------

inline auto GenericVector::rank() const -> size_t
{
  return 1;
}

//-----------------------------------------------------------------------------

inline auto GenericVector::size( size_t dim ) const -> size_t
{
  dolfin_assert( dim == 0 );
  MAYBE_UNUSED( dim );
  return size();
}

//-----------------------------------------------------------------------------

inline void GenericVector::get( real *                 block,
                                const size_t *         num_rows,
                                const size_t * const * rows ) const
{
  get( block, num_rows[0], rows[0] );
}

//-----------------------------------------------------------------------------

inline void GenericVector::set( const real *           block,
                                const size_t *         num_rows,
                                const size_t * const * rows )
{
  set( block, num_rows[0], rows[0] );
}

//-----------------------------------------------------------------------------

inline void GenericVector::add( const real *           block,
                                const size_t *         num_rows,
                                const size_t * const * rows )
{
  add( block, num_rows[0], rows[0] );
}

//-----------------------------------------------------------------------------

inline auto GenericVector::operator[]( size_t i ) const -> real
{
  real value( 0 );
  get( &value, 1, &i );
  return value;
}

//-----------------------------------------------------------------------------

inline auto GenericVector::getitem( size_t i ) const -> real
{
  real value( 0 );
  get( &value, 1, &i );
  return value;
}

//-----------------------------------------------------------------------------

inline void GenericVector::setitem( size_t i, real value )
{
  set( &value, 1, &i );
}

//-----------------------------------------------------------------------------

}

#endif
