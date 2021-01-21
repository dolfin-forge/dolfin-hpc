// Copyright (C) 2006-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_GENERIC_MATRIX_H
#define __DOLFIN_GENERIC_MATRIX_H

#include <dolfin/common/Array.h>
#include <dolfin/la/GenericTensor.h>

namespace dolfin
{

class GenericVector;

/// This class defines a common interface for matrices.

class GenericMatrix : public GenericTensor
{
public:
  /// Destructor
  ~GenericMatrix() override = default;

  //--- Implementation of the GenericTensor interface ---

  /// Initialize zero tensor using sparsity pattern
  void init( const GenericSparsityPattern & sparsity_pattern ) override = 0;

  /// Return copy of tensor
  auto copy() const -> GenericMatrix * override = 0;

  /// Return tensor rank (number of dimensions)
  inline auto rank() const -> uint override;

  /// Return size of given dimension
  auto size( uint dim ) const -> uint override = 0;

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

  //--- Matrix interface ---

  /// Initialize M x N matrix
  virtual void init( uint M, uint N ) = 0;

  /// Initialize M x N matrix
  virtual void init( uint M, uint N, bool distributed ) = 0;

  /// Get block of values
  virtual void get( real * block,
                    uint m, const uint * rows,
                    uint n, const uint * cols ) const = 0;

  /// Set block of values
  virtual void set( const real * block,
                    uint m, const uint * rows,
                    uint n, const uint * cols ) = 0;

  /// Add block of values
  virtual void add( const real * block,
                    uint m, const uint * rows,
                    uint n, const uint * cols ) = 0;

  /// Return norm of matrix
  virtual auto norm( std::string norm_type = "frobenius" ) const -> real = 0;

  /// Get non-zero values of given row
  virtual void getrow( uint            row,
                       Array< uint > & columns,
                       Array< real > & values ) const = 0;

  /// Set values for given row
  virtual void setrow( uint                  row,
                       const Array< uint > & columns,
                       const Array< real > & values ) = 0;

  /// Set given rows to zero
  virtual void zero( uint m, const uint * rows ) = 0;

  /// Set given rows to identity matrix
  virtual void ident( uint m, const uint * rows ) = 0;

  /// Matrix-vector product, y = Ax
  virtual void mult( const GenericVector & x,
                     GenericVector &       y,
                     bool                  transposed = false ) const = 0;

  /// Multiply matrix by given number
  virtual auto operator*=( real a ) -> const GenericMatrix & = 0;

  /// Divide matrix by given number
  virtual auto operator/=( real a ) -> const GenericMatrix & = 0;

  /// Assignment operator
  virtual auto operator=( const GenericMatrix & x ) -> const GenericMatrix & = 0;

  /// Get number of non-zeros in the matrix
  virtual auto nz() const -> uint = 0;

  //--- Convenience functions ---

  /// Get value of given entry
  virtual auto operator()( uint i, uint j ) const -> real;

  /// Get value of given entry
  virtual auto getitem( std::pair< uint, uint > ij ) const -> real;

  /// Set given entry to value
  virtual void setitem( std::pair< uint, uint > ij, real value );
};

//-----------------------------------------------------------------------------
inline auto GenericMatrix::rank() const -> uint
{
  return 2;
}

//-----------------------------------------------------------------------------
inline void GenericMatrix::get( real * block,
                                const uint * num_rows,
                                const uint * const * rows ) const
{
  get( block, num_rows[0], rows[0], num_rows[1], rows[1] );
}

//-----------------------------------------------------------------------------
inline void GenericMatrix::set( const real * block,
                                const uint * num_rows,
                                const uint * const * rows )
{
  set( block, num_rows[0], rows[0], num_rows[1], rows[1] );
}

//-----------------------------------------------------------------------------
inline void GenericMatrix::add( const real * block,
                                const uint * num_rows,
                                const uint * const * rows )
{
  add( block, num_rows[0], rows[0], num_rows[1], rows[1] );
}

//-----------------------------------------------------------------------------
inline auto GenericMatrix::operator()( uint i, uint j ) const -> real
{
  real value( 0 );
  get( &value, 1, &i, 1, &j );
  return value;
}

//-----------------------------------------------------------------------------
inline auto GenericMatrix::getitem( std::pair< uint, uint > ij ) const -> real
{
  real value( 0 );
  get( &value, 1, &ij.first, 1, &ij.second );
  return value;
}

//-----------------------------------------------------------------------------
inline void GenericMatrix::setitem( std::pair< uint, uint > ij, real value )
{
  set( &value, 1, &ij.first, 1, &ij.second );
}

//-----------------------------------------------------------------------------

}

#endif
