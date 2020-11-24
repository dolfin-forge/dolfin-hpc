// Copyright (C) 2006-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_GENERIC_VECTOR_H
#define __DOLFIN_GENERIC_VECTOR_H

#include "GenericSparsityPattern.h"
#include "GenericTensor.h"
#include "VectorNormType.h"
#include "VectorPointwiseOp.h"

#include <dolfin/common/maybe_unused.h>

namespace dolfin
{

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
  GenericVector * copy() const override = 0;

  /// Return tensor rank (number of dimensions)
  uint rank() const override;

  /// Return size of given dimension
  inline uint size( uint dim ) const override;

  /// Get block of values
  inline void get( real * block,
                   const uint * num_rows,
                   const uint * const * rows ) const override;

  /// Set block of values
  inline void set( const real * block,
                   const uint * num_rows,
                   const uint * const * rows ) override;

  /// Add block of values
  inline void add( const real * block,
                   const uint * num_rows,
                   const uint * const * rows ) override;

  /// Set all entries to zero and keep any sparse structure
  void zero() override = 0;

  /// Finalize assembly of tensor
  void apply( FinalizeType finaltype = FINALIZE ) override = 0;

  /// Display tensor
  void disp( uint precision = 2 ) const override = 0;

  //--- Vector interface ---

  /// Initialize vector of size N
  virtual void init( uint N ) = 0;

  /// Initialize vector of size N and distribute if specified
  virtual void init( uint N, bool distributed ) = 0;

  /// Initialize ghost entries
  virtual void init_ghosted( uint                         n,
                             _ordered_set< uint > &       indices,
                             _ordered_map< uint, uint > & map ) = 0;

  /// Return size of vector
  virtual uint size() const = 0;

  /// Return local size of vector
  virtual uint local_size() const = 0;

  /// Return rank's offset into vector
  virtual uint offset() const = 0;

  /// Get block of values
  virtual void get( real * block, uint m, const uint * rows ) const = 0;

  /// Set block of values
  virtual void set( const real * block, uint m, const uint * rows ) = 0;

  /// Add block of values
  virtual void add( const real * block, uint m, const uint * rows ) = 0;

  /// Get all values
  virtual void get( real * values ) const = 0;

  /// Set all values
  virtual void set( real * values ) = 0;

  /// Add values to each entry
  virtual void add( real * values ) = 0;

  /// Add multiple of given vector (y=a*x+y)
  virtual void axpy( real a, const GenericVector & x ) = 0;

  /// Add multiple of given vector (y=a*x+b*y)
  virtual void axpby( real a, const GenericVector & x,
                      real b ) = 0;

  /// Add multiple of given vector (w=a*x+y)
  virtual void waxpy( real a, const GenericVector & x,
                              const GenericVector & y ) = 0;

  /// Add multiple of given vector (z=a*x+b*y+c*z)
  virtual void axpbypcz( real a, const GenericVector & x,
                         real b, const GenericVector & y,
                         real c ) = 0;

  /// Return inner product with given vector
  virtual real inner( const GenericVector & x ) const = 0;

  /// Return norm of vector
  virtual real norm( VectorNormType type = l2 ) const = 0;

  /// Return minimum value of vector
  virtual real min() const = 0;

  /// Return maximum value of vector
  virtual real max() const = 0;

  /// Return pointwise operator op of vector and given vector x
  virtual void pointwise( const GenericVector & x,
                          VectorPointwiseOp     op = pw_min ) const = 0;

  /// Multiply vector by given number
  virtual const GenericVector & operator*=( real a ) = 0;

  /// Divide vector by given number
  virtual const GenericVector & operator/=( real a ) = 0;

  /// Multiply vector by given vector component-wise
  virtual const GenericVector & operator*=( const GenericVector & x ) = 0;

  /// Add given vector
  virtual const GenericVector & operator+=( const GenericVector & x ) = 0;

  /// Subtract given vector
  virtual const GenericVector & operator-=( const GenericVector & x ) = 0;

  /// Assignment operator
  virtual const GenericVector & operator=( const GenericVector & x ) = 0;

  /// Assignment operator
  virtual const GenericVector & operator=( real a ) = 0;

  //--- Convenience functions ---

  /// Get value of given entry
  virtual real operator[]( uint i ) const;

  /// Get value of given entry
  virtual real getitem( uint i ) const;

  /// Set given entry to value
  virtual void setitem( uint i, real value );
};

//-----------------------------------------------------------------------------
inline void GenericVector::init( const GenericSparsityPattern & sparsity_pattern )
{
  init( sparsity_pattern.size( 0 ), sparsity_pattern.is_distributed() );
}

//-----------------------------------------------------------------------------
inline uint GenericVector::rank() const
{
  return 1;
}

//-----------------------------------------------------------------------------
inline uint GenericVector::size( uint dim ) const
{
  dolfin_assert( dim == 0 );
  MAYBE_UNUSED( dim );
  return size();
}

//-----------------------------------------------------------------------------
inline void GenericVector::get( real *               block,
                                const uint *         num_rows,
                                const uint * const * rows ) const
{
  get( block, num_rows[0], rows[0] );
}

//-----------------------------------------------------------------------------
inline void GenericVector::set( const real *         block,
                                const uint *         num_rows,
                                const uint * const * rows )
{
  set( block, num_rows[0], rows[0] );
}

//-----------------------------------------------------------------------------
inline void GenericVector::add( const real *         block,
                                const uint *         num_rows,
                                const uint * const * rows )
{
  add( block, num_rows[0], rows[0] );
}

//-----------------------------------------------------------------------------
inline real GenericVector::operator[]( uint i ) const
{
  real value( 0 );
  get( &value, 1, &i );
  return value;
}

//-----------------------------------------------------------------------------
inline real GenericVector::getitem( uint i ) const
{
  real value( 0 );
  get( &value, 1, &i );
  return value;
}

//-----------------------------------------------------------------------------
inline void GenericVector::setitem( uint i, real value )
{
  set( &value, 1, &i );
}

//-----------------------------------------------------------------------------

}

#endif
