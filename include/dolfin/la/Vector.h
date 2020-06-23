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
  Vector()
    : vector_( DefaultFactory::factory().createVector() )
  {
  }

  /// Create vector of size N distributed by default
  explicit Vector( uint N )
    : vector_( DefaultFactory::factory().createVector() )
  {
    vector_->init( N );
  }

  /// Create vector of size N distributed if specified
  explicit Vector( uint N, bool distributed )
    : vector_( DefaultFactory::factory().createVector() )
  {
    vector_->init( N, distributed );
  }

  /// Copy constructor
  explicit Vector( Vector const & x )
    : vector_( x.vector_->copy() )
  {
  }

  /// Destructor
  ~Vector()
  {
    delete vector_;
  }

  //--- Implementation of the GenericTensor interface ---

  /// Return copy of tensor
  Vector * copy() const;

  /// Set all entries to zero and keep any sparse structure
  void zero();

  /// Finalize assembly of tensor
  void apply( FinalizeType finaltype = FINALIZE );

  /// Display tensor
  void disp( uint precision = 2 ) const;

  //--- Implementation of the GenericVector interface ---

  /// Initialize vector of size N
  void init( uint N );

  /// Initialize vector of size N and distribute if specified
  void init( uint N, bool distributed );

  void init_ghosted( uint                         n,
                     _ordered_set< uint > &       indices,
                     _ordered_map< uint, uint > & map );

  /// Return size of vector
  uint size() const;

  /// Return local size of vector
  uint local_size() const;

  /// Return rank's offset into vector
  uint offset() const;

  /// Get block of values
  void get( real * block, uint m, const uint * rows ) const;

  /// Set block of values
  void set( const real * block, uint m, const uint * rows );

  /// Add block of values
  void add( const real * block, uint m, const uint * rows );

  /// Get all local values (not ghost entries)
  void get( real * values ) const;

  /// Set all local values (not ghost entries)
  void set( real * values );

  /// Add values to each local entry (not ghost entries)
  void add( real * values );

  /// Add multiple of given vector (AXPY operation)
  void axpy( real a, const GenericVector & x );

  /// Return inner product with given vector
  real inner( const GenericVector & x ) const;

  /// Return norm of vector
  real norm( VectorNormType type = l2 ) const;

  /// Return minimum value of vector
  real min() const;

  /// Return maximum value of vector
  real max() const;

  /// Return pointwise operator op of vector and given vector x
  void pointwise( const GenericVector & x,
                  VectorPointwiseOp     op = pw_min ) const;

  /// Multiply vector by given number
  Vector & operator*=( real a );

  /// Divide vector by given number
  Vector & operator/=( real a );

  /// Multiply vector by given vector component-wise
  Vector & operator*=( const GenericVector & x );

  /// Add given vector
  Vector & operator+=( const GenericVector & x );

  /// Subtract given vector
  Vector & operator-=( const GenericVector & x );

  /// Assignment operator
  GenericVector & operator=( const GenericVector & x );

  /// Assignment operator
  Vector & operator=( real a );

  //--- Special functions ---

  /// Return linear algebra backend factory
  LinearAlgebraFactory & factory() const;

  //--- Special functions, intended for library use only ---

  /// Return concrete instance / unwrap (const)
  GenericVector const * instance() const;

  /// Return concrete instance / unwrap (non-const version)
  GenericVector * instance();

  //--- Special Vector functions ---

  /// Assignment operator
  Vector & operator=( Vector const & x );

private:
  // Pointer to concrete implementation
  GenericVector * const vector_;
};

//-----------------------------------------------------------------------------
inline Vector * Vector::copy() const
{
  return new Vector( *this );
}

//-----------------------------------------------------------------------------
inline void Vector::zero()
{
  vector_->zero();
}

//-----------------------------------------------------------------------------
inline void Vector::apply( FinalizeType finaltype )
{
  vector_->apply( finaltype );
}

//-----------------------------------------------------------------------------
inline void Vector::disp( uint precision ) const
{
  vector_->disp( precision );
}

//-----------------------------------------------------------------------------
inline void Vector::init( uint N )
{
  vector_->init( N );
}

//-----------------------------------------------------------------------------
inline void Vector::init( uint N, bool distributed )
{
  vector_->init( N, distributed );
}

//-----------------------------------------------------------------------------
inline void Vector::init_ghosted( uint                         n,
                                  _ordered_set< uint > &       indices,
                                  _ordered_map< uint, uint > & map )
{
  vector_->init_ghosted( n, indices, map );
}

//-----------------------------------------------------------------------------
inline uint Vector::size() const
{
  return vector_->size();
}

//-----------------------------------------------------------------------------
inline uint Vector::local_size() const
{
  return vector_->local_size();
}

//-----------------------------------------------------------------------------
inline uint Vector::offset() const
{
  return vector_->offset();
}

//-----------------------------------------------------------------------------
inline void Vector::get( real * block, uint m, const uint * rows ) const
{
  vector_->get( block, m, rows );
}

//-----------------------------------------------------------------------------
inline void Vector::set( const real * block, uint m, const uint * rows )
{
  vector_->set( block, m, rows );
}

//-----------------------------------------------------------------------------
inline void Vector::add( const real * block, uint m, const uint * rows )
{
  vector_->add( block, m, rows );
}

//-----------------------------------------------------------------------------
inline void Vector::get( real * values ) const
{
  vector_->get( values );
}

//-----------------------------------------------------------------------------
inline void Vector::set( real * values )
{
  vector_->set( values );
}

//-----------------------------------------------------------------------------
inline void Vector::add( real * values )
{
  vector_->add( values );
}

//-----------------------------------------------------------------------------
inline void Vector::axpy( real a, const GenericVector & x )
{
  vector_->axpy( a, x );
}

//-----------------------------------------------------------------------------
inline real Vector::inner( const GenericVector & x ) const
{
  return vector_->inner( x );
}

//-----------------------------------------------------------------------------
inline real Vector::norm( VectorNormType type ) const
{
  return vector_->norm( type );
}

//-----------------------------------------------------------------------------
inline real Vector::min() const
{
  return vector_->min();
}

//-----------------------------------------------------------------------------
inline real Vector::max() const
{
  return vector_->max();
}

//-----------------------------------------------------------------------------
inline void Vector::pointwise( const GenericVector & x,
                               VectorPointwiseOp     op ) const
{
  return vector_->pointwise( x, op );
}

//-----------------------------------------------------------------------------
inline Vector & Vector::operator*=( real a )
{
  *vector_ *= a;
  return *this;
}

//-----------------------------------------------------------------------------
inline Vector & Vector::operator/=( real a )
{
  *this *= 1.0 / a;
  return *this;
}

//-----------------------------------------------------------------------------
inline Vector & Vector::operator*=( const GenericVector & x )
{
  *vector_ *= x;
  return *this;
}

//-----------------------------------------------------------------------------
inline Vector & Vector::operator+=( const GenericVector & x )
{
  axpy( 1.0, x );
  return *this;
}

//-----------------------------------------------------------------------------
inline Vector & Vector::operator-=( const GenericVector & x )
{
  axpy( -1.0, x );
  return *this;
}

//-----------------------------------------------------------------------------
inline GenericVector & Vector::operator=( const GenericVector & x )
{
  *vector_ = x;
  return *this;
}

//-----------------------------------------------------------------------------
inline Vector & Vector::operator=( real a )
{
  *vector_ = a;
  return *this;
}

//-----------------------------------------------------------------------------
inline LinearAlgebraFactory & Vector::factory() const
{
  return vector_->factory();
}

//-----------------------------------------------------------------------------
inline GenericVector const * Vector::instance() const
{
  return vector_;
}

//-----------------------------------------------------------------------------
inline GenericVector * Vector::instance()
{
  return vector_;
}

//-----------------------------------------------------------------------------
inline Vector & Vector::operator=( Vector const & x )
{
  *vector_ = *x.vector_;
  return *this;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_VECTOR_H */
