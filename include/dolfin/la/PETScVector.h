// Copyright (C) 2004-2008 Johan Hoffman, Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PETSC_VECTOR_H
#define __DOLFIN_PETSC_VECTOR_H

#include <dolfin/common/Variable.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/log/LogStream.h>

#ifdef HAVE_PETSC

#include <dolfin/la/PETScObject.h>

#include <petscvec.h>

namespace dolfin
{

/// This class provides a simple vector class based on PETSc.
/// It is a simple wrapper for a PETSc vector pointer (Vec)
/// implementing the GenericVector interface.
///
/// The interface is intentionally simple. For advanced usage,
/// access the PETSc Vec pointer using the function vec() and
/// use the standard PETSc interface.

class PETScVector : public GenericVector, public PETScObject, public Variable
{

public:
  /// Create empty vector
  PETScVector();

  /// Create vector of local size N
  explicit PETScVector( size_t N, bool distributed = true );

  /// Copy constructor
  explicit PETScVector( const PETScVector & x );

  /// Create vector from given PETSc Vec pointer
  explicit PETScVector( Vec x );

  /// Destructor
  ~PETScVector() override;

  //--- Implementation of the GenericTensor interface ---

  /// Return copy of tensor
  auto copy() const -> PETScVector * override;

  /// Set all entries to zero and keep any sparse structure
  void zero() override;

  /// Finalize assembly of tensor
  void apply( FinalizeType finaltype = FINALIZE ) override;

  /// Display tensor
  void disp( size_t precision = 0 ) const override;

  //--- Implementation of the GenericVector interface ---

  /// Initialize vector of local size N, distributed by default
  void init( size_t N ) override;

  /// Initialize vector of local size N, distributed if specified
  void init( size_t N, bool distributed ) override;

  /// Initialize ghost entries
  void init_ghosted( size_t                           n,
                     _ordered_set< size_t > &         indices,
                     _ordered_map< size_t, size_t > & map ) override;

  /// Return size of vector
  auto size() const -> size_t override;

  /// Return local size of vector
  auto local_size() const -> size_t override;

  /// Return rank's offset into vector
  auto offset() const -> size_t override;

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
  auto inner( const GenericVector & v ) const -> real override;

  /// Return norm of vector
  auto norm( VectorNormType type = l2 ) const -> real override;

  /// Return minimum value of vector
  auto min() const -> real override;

  /// Return maximum value of vector
  auto max() const -> real override;

  /// Return pointwise operator op of vector and given vector x
  void pointwise( const GenericVector & x,
                  VectorPointwiseOp     op = pw_min ) const override;

  /// Multiply vector by given number
  auto operator*=( real a ) -> PETScVector & override;

  /// Divide vector by given number
  auto operator/=( real a ) -> PETScVector & override;

  /// Multiply vector by given vector component-wise
  auto operator*=( const GenericVector & x ) -> PETScVector & override;

  /// Add given vector
  auto operator+=( const GenericVector & x ) -> PETScVector & override;

  /// Subtract given vector
  auto operator-=( const GenericVector & x ) -> PETScVector & override;

  /// Assignment operator
  auto operator=( const GenericVector & x ) -> PETScVector & override;

  /// Assignment operator
  auto operator=( real a ) -> PETScVector & override;

  //--- Special functions ---

  /// Return linear algebra backend factory
  auto factory() const -> LinearAlgebraFactory & override;

  //--- Special PETSc functions ---

  /// Assignment operator
  auto operator=( PETScVector const & x ) -> PETScVector &;

  /// Return PETSc Vec pointer
  auto vec() const -> Vec;

  inline auto ghosted() -> bool
  {
    return is_ghosted_;
  }

private:
  //
  void clear();

  // PETSc Vec pointer
  Vec x_ { nullptr };

  // True if the vector is distributed
  bool is_distributed_ { false };

  // True if the vector has ghost points
  bool is_ghosted_ { false };

  using GhostMapping = _map< int, int >;
  GhostMapping mapping_;
};

//-----------------------------------------------------------------------------
inline void PETScVector::clear()
{
  if ( x_ )
  {
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    VecDestroy( &x_ );
#else
    VecDestroy( x_ );
#endif
    is_ghosted_     = false;
    is_distributed_ = false;
  }
}
//-----------------------------------------------------------------------------
inline void PETScVector::init( size_t N )
{
  init( N, true );
}
//-----------------------------------------------------------------------------
inline auto PETScVector::copy() const -> PETScVector *
{
  return new PETScVector( *this );
}
//-----------------------------------------------------------------------------
inline void
  PETScVector::set( const real * block, size_t m, const size_t * rows )
{
  dolfin_assert( x_ );

  PetscInt M_ = m;

  // FIXME this is potentially costly
  std::vector< PetscInt > rows_( rows, rows + m );

  VecSetValues( x_, M_, rows_.data(), block, INSERT_VALUES );
}
//-----------------------------------------------------------------------------
inline void
  PETScVector::add( const real * block, size_t m, const size_t * rows )
{
  dolfin_assert( x_ );

  PetscInt M_ = m;

  // FIXME this is potentially costly
  std::vector< PetscInt > rows_( rows, rows + m );

  VecSetValues( x_, M_, rows_.data(), block, ADD_VALUES );
}
//-----------------------------------------------------------------------------
inline void PETScVector::apply( FinalizeType )
{

  VecAssemblyBegin( x_ );
  VecAssemblyEnd( x_ );

  if ( is_ghosted_ )
  {
    VecGhostUpdateBegin( x_, INSERT_VALUES, SCATTER_FORWARD );
    VecGhostUpdateEnd( x_, INSERT_VALUES, SCATTER_FORWARD );
  }
}
//-----------------------------------------------------------------------------
inline void PETScVector::zero()
{
  dolfin_assert( x_ );
  real a = 0.0;
  VecSet( x_, a );
}
//-----------------------------------------------------------------------------
inline auto PETScVector::size() const -> size_t
{
  if ( x_ == nullptr )
    return 0;
  PetscInt n = 0;
  VecGetSize( x_, &n );
  return static_cast< size_t >( n );
}
//-----------------------------------------------------------------------------
inline auto PETScVector::local_size() const -> size_t
{
  dolfin_assert( x_ );
  PetscInt n = 0;
  VecGetLocalSize( x_, &n );
  return static_cast< size_t >( n );
}
//-----------------------------------------------------------------------------
inline auto PETScVector::offset() const -> size_t
{
  dolfin_assert( x_ );
  PetscInt low, high;
  VecGetOwnershipRange( x_, &low, &high );
  return static_cast< size_t >( low );
}
//-----------------------------------------------------------------------------
inline auto PETScVector::operator=( const GenericVector & v ) -> PETScVector &
{
  *this = v.down_cast< PETScVector >();
  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScVector::operator=( PETScVector const & v ) -> PETScVector &
{
  if ( &v != this )
  {
    dolfin_assert( v.x_ );
    init( v.local_size(), v.is_distributed_ );
    VecCopy( v.x_, x_ );
  }
  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScVector::operator=( real a ) -> PETScVector &
{
  dolfin_assert( x_ );
  VecSet( x_, a );
  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScVector::operator*=( const GenericVector & y ) -> PETScVector &
{
  dolfin_assert( x_ );
  PETScVector const & v = y.down_cast< PETScVector >();
  dolfin_assert( v.x_ );

  if ( size() != v.size() )
  {
    error(
      "Vectors must have the same size for componentwise multiplication." );
  }

  VecPointwiseMult( x_, x_, v.x_ );

  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScVector::operator+=( const GenericVector & x ) -> PETScVector &
{
  this->axpy( 1.0, x );
  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScVector::operator-=( const GenericVector & x ) -> PETScVector &
{
  this->axpy( -1.0, x );
  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScVector::operator*=( const real a ) -> PETScVector &
{
  dolfin_assert( x_ );
  VecScale( x_, a );

  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScVector::operator/=( const real a ) -> PETScVector &
{
  dolfin_assert( x_ );
  dolfin_assert( a != 0.0 );

  const real b = 1.0 / a;
  VecScale( x_, b );

  return *this;
}
//-----------------------------------------------------------------------------
inline auto PETScVector::inner( const GenericVector & y ) const -> real
{
  dolfin_assert( x_ );

  PETScVector const & v = y.down_cast< PETScVector >();
  dolfin_assert( v.x_ );

  real a;
  VecDot( v.x_, x_, &a );

  return a;
}
//-----------------------------------------------------------------------------
inline void PETScVector::axpy( real a, const GenericVector & x )
{
  dolfin_assert( x_ );

  PETScVector const & v = x.down_cast< PETScVector >();
  dolfin_assert( v.x_ );

  if ( size() != v.size() )
  {
    error( "The vectors must be of the same size to apply AXPY." );
  }

  VecAXPY( x_, a, v.x_ );
}
//-----------------------------------------------------------------------------
inline void PETScVector::axpby( real a, const GenericVector & x, real b )
{
  dolfin_assert( x_ );

  PETScVector const & v = x.down_cast< PETScVector >();
  dolfin_assert( v.x_ );

  if ( size() != v.size() )
  {
    error( "The vectors must be of the same size to apply AXPBY." );
  }

  VecAXPBY( x_, a, b, v.x_ );
}
//-----------------------------------------------------------------------------
inline void
  PETScVector::waxpy( real a, const GenericVector & x, const GenericVector & y )
{
  dolfin_assert( x_ );

  PETScVector const & v = x.down_cast< PETScVector >();
  dolfin_assert( v.x_ );

  PETScVector const & u = y.down_cast< PETScVector >();
  dolfin_assert( u.x_ );

  if ( size() != v.size() or size() != u.size() )
  {
    error( "The vectors must be of the same size to apply WAXPY." );
  }

  VecWAXPY( x_, a, v.x_, u.x_ );
}
//-----------------------------------------------------------------------------
inline void PETScVector::axpbypcz( real                  a,
                                   const GenericVector & x,
                                   real                  b,
                                   const GenericVector & y,
                                   real                  c )
{
  dolfin_assert( x_ );

  PETScVector const & v = x.down_cast< PETScVector >();
  dolfin_assert( v.x_ );

  PETScVector const & u = y.down_cast< PETScVector >();
  dolfin_assert( u.x_ );

  if ( size() != v.size() or size() != u.size() )
  {
    error( "The vectors must be of the same size to apply VecAXPBYPCZ." );
  }

  VecAXPBYPCZ( x_, a, b, c, v.x_, u.x_ );
}
//-----------------------------------------------------------------------------
inline auto PETScVector::min() const -> real
{
  real value = 0.0;

  VecMin( x_, PETSC_NULL, &value );

  return value;
}
//-----------------------------------------------------------------------------
inline auto PETScVector::max() const -> real
{
  real value = 0.0;

  VecMax( x_, PETSC_NULL, &value );

  return value;
}

//-----------------------------------------------------------------------------
inline auto PETScVector::vec() const -> Vec
{
  return x_;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* HAVE_PETSC */

#endif /* __DOLFIN_PETSC_VECTOR_H */
