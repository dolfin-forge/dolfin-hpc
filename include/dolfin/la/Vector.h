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
  explicit Vector( size_t N )
    : vector_( DefaultFactory::factory().createVector() )
  {
    vector_->init( N );
  }

  /// Create vector of size N distributed if specified
  explicit Vector( size_t N, bool distributed )
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
  ~Vector() override
  {
    delete vector_;
  }

  //--- Implementation of the GenericTensor interface ---

  /// Return copy of tensor
  Vector * copy() const override;

  /// Set all entries to zero and keep any sparse structure
  void zero() override;

  /// Finalize assembly of tensor
  void apply( FinalizeType finaltype = FINALIZE ) override;

  /// Display tensor
  void disp( size_t precision = 2 ) const override;

  //--- Implementation of the GenericVector interface ---

  void init( const GenericSparsityPattern & sparsity_pattern ) override;

  /// Initialize vector of size N
  void init( size_t N ) override;

  /// Initialize vector of size N and distribute if specified
  void init( size_t N, bool distributed ) override;

  void init_ghosted( size_t                           n,
                     _ordered_set< size_t > &         indices,
                     _ordered_map< size_t, size_t > & map ) override;

  /// Return size of vector
  size_t size() const override;

  /// Return local size of vector
  size_t local_size() const override;

  /// Return rank's offset into vector
  size_t offset() const override;

  /// Get block of values
  void get( real * block, size_t m, const size_t * rows ) const override;

  /// Set block of values
  void set( const real * block, size_t m, const size_t * rows ) override;

  /// Add block of values
  void add( const real * block, size_t m, const size_t * rows ) override;

  /// Get all local values (not ghost entries)
  void get( real * values ) const override;

  /// Set all local values (not ghost entries)
  void set( real * values ) override;

  /// Add values to each local entry (not ghost entries)
  void add( real * values ) override;

  /// Add multiple of given vector (y=a*x+y)
  void axpy( real a, const GenericVector & x ) override;

  /// Add multiple of given vector (y=a*x+b*y)
  void axpby( real a, const GenericVector & x, real b ) override;

  /// Add multiple of given vector (w=a*x+y)
  void
    waxpy( real a, const GenericVector & x, const GenericVector & y ) override;

  /// Add multiple of given vector (z=a*x+b*y+c*z)
  void axpbypcz( real                  a,
                 const GenericVector & x,
                 real                  b,
                 const GenericVector & y,
                 real                  c ) override;

  /// Return inner product with given vector
  real inner( const GenericVector & x ) const override;

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
  GenericVector & operator=( const GenericVector & x ) override;

  /// Assignment operator
  Vector & operator=( real a ) override;

  //--- Special functions ---

  /// Return linear algebra backend factory
  LinearAlgebraFactory & factory() const override;

  //--- Special functions, intended for library use only ---

  /// Return concrete instance / unwrap (const)
  GenericVector const * instance() const override;

  /// Return concrete instance / unwrap (non-const version)
  GenericVector * instance() override;

  //--- Special Vector functions ---

  /// Assignment operator
  Vector & operator=( Vector const & x );

private:
  // Pointer to concrete implementation
  GenericVector * vector_ { nullptr };
};

//-----------------------------------------------------------------------------

inline Vector * Vector::copy() const
{
  Vector * V = new Vector();
  delete V->vector_;
  V->vector_ = vector_->copy();
  return V;
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

inline void Vector::disp( size_t precision ) const
{
  vector_->disp( precision );
}

//-----------------------------------------------------------------------------


inline void Vector::init( const GenericSparsityPattern & sparsity_pattern )
{
  vector_->init( sparsity_pattern );
}

//-----------------------------------------------------------------------------

inline void Vector::init( size_t N )
{
  vector_->init( N );
}

//-----------------------------------------------------------------------------

inline void Vector::init( size_t N, bool distributed )
{
  vector_->init( N, distributed );
}

//-----------------------------------------------------------------------------

inline void Vector::init_ghosted( size_t                           n,
                                  _ordered_set< size_t > &         indices,
                                  _ordered_map< size_t, size_t > & map )
{
  vector_->init_ghosted( n, indices, map );
}

//-----------------------------------------------------------------------------

inline size_t Vector::size() const
{
  return vector_->size();
}

//-----------------------------------------------------------------------------

inline size_t Vector::local_size() const
{
  return vector_->local_size();
}

//-----------------------------------------------------------------------------

inline size_t Vector::offset() const
{
  return vector_->offset();
}

//-----------------------------------------------------------------------------

inline void Vector::get( real * block, size_t m, const size_t * rows ) const
{
  vector_->get( block, m, rows );
}

//-----------------------------------------------------------------------------

inline void Vector::set( const real * block, size_t m, const size_t * rows )
{
  vector_->set( block, m, rows );
}

//-----------------------------------------------------------------------------

inline void Vector::add( const real * block, size_t m, const size_t * rows )
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

inline void Vector::axpby( real a, const GenericVector & x, real b )
{
  vector_->axpby( a, x, b );
}

//-----------------------------------------------------------------------------

inline void
  Vector::waxpy( real a, const GenericVector & x, const GenericVector & y )
{
  vector_->waxpy( a, x, y );
}

//-----------------------------------------------------------------------------

inline void Vector::axpbypcz( real                  a,
                              const GenericVector & x,
                              real                  b,
                              const GenericVector & y,
                              real                  c )
{
  vector_->axpbypcz( a, x, b, y, c );
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
