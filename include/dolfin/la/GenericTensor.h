// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_GENERIC_TENSOR_H
#define __DOLFIN_GENERIC_TENSOR_H

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

enum FinalizeType
{
  FINALIZE,
  FLUSH
};

namespace dolfin
{

class GenericSparsityPattern;
class LinearAlgebraFactory;

//-----------------------------------------------------------------------------

/// This class defines a common interface for arbitrary rank tensors.

class GenericTensor
{
public:
  /// Destructor
  virtual ~GenericTensor() = default;

  //--- Basic GenericTensor interface ---

  /// Initialize zero tensor using sparsity pattern
  virtual void init( const GenericSparsityPattern & sparsity_pattern ) = 0;

  /// Return copy of tensor
  virtual auto copy() const -> GenericTensor * = 0;

  /// Return tensor rank (number of dimensions)
  virtual auto rank() const -> size_t = 0;

  /// Return size of given dimension
  virtual auto size( size_t dim ) const -> size_t = 0;

  /// Get block of values
  virtual void get( real *                 block,
                    const size_t *         num_rows,
                    const size_t * const * rows ) const = 0;

  /// Set block of values
  virtual void set( const real *           block,
                    const size_t *         num_rows,
                    const size_t * const * rows ) = 0;

  /// Add block of values
  virtual void add( const real *           block,
                    const size_t *         num_rows,
                    const size_t * const * rows ) = 0;

  /// Set all entries to zero and keep any sparse structure
  virtual void zero() = 0;

  /// Finalize assembly of tensor
  virtual void apply( FinalizeType finaltype = FINALIZE ) = 0;

  /// Display tensor
  virtual void disp( size_t precision = 2 ) const = 0;

  //--- Special functions, downcasting to concrete types ---

  /// Return linear algebra backend factory
  virtual auto factory() const -> LinearAlgebraFactory & = 0;

  /// Cast a GenericTensor to its derived class (const)
  template < class T >
  auto down_cast() const -> const T &;

  /// Cast a GenericTensor to its derived class (non-const version)
  template < class T >
  auto down_cast() -> T &;

  /// Check whether the GenericTensor instance matches a specific type
  template < class T >
  auto has_type() const -> bool;

  //--- Special functions, intended for library use only ---

  /// Return concrete instance / unwrap (const)
  virtual auto instance() const -> const GenericTensor *;

  /// Return concrete instance / unwrap (non-const version)
  virtual auto instance() -> GenericTensor *;
};


//-----------------------------------------------------------------------------

template < class T >
inline auto GenericTensor::down_cast() const -> const T &
{
  const T * t = dynamic_cast< const T * >( instance() );
  if ( !t )
    error( "GenericTensor cannot be cast to the requested type." );
  return *t;
}

//-----------------------------------------------------------------------------

template < class T >
inline auto GenericTensor::down_cast() -> T &
{
  T * t = dynamic_cast< T * >( instance() );
  if ( !t )
    error( "GenericTensor cannot be cast to the requested type." );
  return *t;
}

//-----------------------------------------------------------------------------

template < class T >
inline auto GenericTensor::has_type() const -> bool
{
  return bool( dynamic_cast< const T * >( instance() ) );
}

//-----------------------------------------------------------------------------

inline auto GenericTensor::instance() const -> const GenericTensor *
{
  return this;
}

//-----------------------------------------------------------------------------

inline auto GenericTensor::instance() -> GenericTensor *
{
  return this;
}

//-----------------------------------------------------------------------------

}

#endif
