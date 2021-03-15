// Copyright (C) 2006-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_GENERIC_MATRIX_H
#define __DOLFIN_GENERIC_MATRIX_H

#include <dolfin/la/GenericTensor.h>

#include <vector>

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
  inline auto rank() const -> size_t override;

  /// Return size of given dimension
  auto size( size_t dim ) const -> size_t override = 0;

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

  //--- Matrix interface ---

  /// Initialize M x N matrix
  virtual void init( size_t M, size_t N ) = 0;

  /// Initialize M x N matrix
  virtual void init( size_t M, size_t N, bool distributed ) = 0;

  /// Get block of values
  virtual void get( real *         block,
                    size_t         m,
                    const size_t * rows,
                    size_t         n,
                    const size_t * cols ) const = 0;

  /// Set block of values
  virtual void set( const real *   block,
                    size_t         m,
                    const size_t * rows,
                    size_t         n,
                    const size_t * cols ) = 0;

  /// Add block of values
  virtual void add( const real *   block,
                    size_t         m,
                    const size_t * rows,
                    size_t         n,
                    const size_t * cols ) = 0;

  /// Return norm of matrix
  virtual auto norm( std::string norm_type = "frobenius" ) const -> real = 0;

  /// Get non-zero values of given row
  virtual void getrow( size_t                  row,
                       std::vector< size_t > & columns,
                       std::vector< real > &   values ) const = 0;

  /// Set values for given row
  virtual void setrow( size_t                        row,
                       const std::vector< size_t > & columns,
                       const std::vector< real > &   values ) = 0;

  /// Set given rows to zero
  virtual void zero( size_t m, const size_t * rows ) = 0;

  /// Set given rows to identity matrix
  virtual void ident( size_t m, const size_t * rows ) = 0;

  /// Matrix-vector product, y = Ax
  virtual void mult( const GenericVector & x,
                     GenericVector &       y,
                     bool                  transposed = false ) const = 0;

  /// Multiply matrix by given number
  virtual auto operator*=( real a ) -> const GenericMatrix & = 0;

  /// Divide matrix by given number
  virtual auto operator/=( real a ) -> const GenericMatrix & = 0;

  /// Assignment operator
  virtual auto operator      =( const GenericMatrix & x )
    -> const GenericMatrix & = 0;

  /// Get number of non-zeros in the matrix
  virtual auto nz() const -> size_t = 0;

  //--- Convenience functions ---

  /// Get value of given entry
  virtual auto operator()( size_t i, size_t j ) const -> real;

  /// Get value of given entry
  virtual auto getitem( std::pair< size_t, size_t > ij ) const -> real;

  /// Set given entry to value
  virtual void setitem( std::pair< size_t, size_t > ij, real value );
};

//-----------------------------------------------------------------------------
inline auto GenericMatrix::rank() const -> size_t
{
  return 2;
}

//-----------------------------------------------------------------------------
inline void GenericMatrix::get( real *                 block,
                                const size_t *         num_rows,
                                const size_t * const * rows ) const
{
  get( block, num_rows[0], rows[0], num_rows[1], rows[1] );
}

//-----------------------------------------------------------------------------
inline void GenericMatrix::set( const real *           block,
                                const size_t *         num_rows,
                                const size_t * const * rows )
{
  set( block, num_rows[0], rows[0], num_rows[1], rows[1] );
}

//-----------------------------------------------------------------------------
inline void GenericMatrix::add( const real *           block,
                                const size_t *         num_rows,
                                const size_t * const * rows )
{
  add( block, num_rows[0], rows[0], num_rows[1], rows[1] );
}

//-----------------------------------------------------------------------------
inline auto GenericMatrix::operator()( size_t i, size_t j ) const -> real
{
  real value( 0 );
  get( &value, 1, &i, 1, &j );
  return value;
}

//-----------------------------------------------------------------------------
inline auto GenericMatrix::getitem( std::pair< size_t, size_t > ij ) const
  -> real
{
  real value( 0 );
  get( &value, 1, &ij.first, 1, &ij.second );
  return value;
}

//-----------------------------------------------------------------------------
inline void GenericMatrix::setitem( std::pair< size_t, size_t > ij, real value )
{
  set( &value, 1, &ij.first, 1, &ij.second );
}

//-----------------------------------------------------------------------------

}

#endif
