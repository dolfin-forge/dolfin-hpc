// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

// This class is now reimplemented in terms of a template.
// Additionally it specifies a MAX_SIZE public attribute to avoid the madness of
// having half of the classes with 2 or 3 hard-coded as dimension.

#ifndef __DOLFIN_POINT_H
#define __DOLFIN_POINT_H

#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Space.h>

#include <iostream>

namespace dolfin
{

/// A point represents a point in R^d with coordinates x, y, z, or,
/// alternatively, a vector in R^d, supporting standard operations
/// like the norm, distances, scalar and vector products etc.

template < uint D = Space::MAX_DIMENSION >
class point
{
public:
  static uint const MAX_SIZE = Space::MAX_DIMENSION;

  /// Create a point at (x, y, z)
  point( const real x = 0.0, const real y = 0.0, const real z = 0.0 )
  {
    x_[0] = x;

    if ( D > 1 )
      x_[1] = y;

    if ( D > 2 )
      x_[2] = z;
  }

  /// Create a point at given coordinates
  point( real const * x )
  {
    std::copy( x, x + D, x_ );
  }

  /// Create a point at given coordinates in R^d
  point( uint d, real const * x )
  {
    std::fill( x_, x_ + D, 0.0 );
    std::copy( x, x + d, x_ );
  }

  /// Copy constructor
  point( point const & p )
  {
    std::copy( p.x_, p.x_ + D, x_ );
  }

  /// Assignment operator
  point & operator=( point const & p )
  {
    std::copy( p.x_, p.x_ + D, x_ );
    return *this;
  }

  /// Destructor
  ~point()
  {
  }

  /// Set coordinates for given Euclidean space dimension
  void set( real const * x, uint dim );

  /// Return coordinate in direction i
  real & operator[]( uint i );

  /// Return coordinate in direction i
  real const & operator[]( uint i ) const;

  /// cast operator
  operator real *();

  /// cast operator
  operator real const *() const;

  /// Addition assignment
  point & operator+=( point const & p );

  /// Addition
  point operator+( point const & p ) const;

  /// Subtraction assignment
  point & operator-=( point const & p );

  /// Subtraction
  point operator-( point const & p ) const;

  /// Scalar multiplication assignment
  point & operator*=( real a );

  /// Scalar multiplication
  point operator*( real a ) const;

  /// Scalar division assignment
  point & operator/=( real a );

  /// Scalar division
  point operator/( real a ) const;

  /// Compute distance to given point
  real dist( point const & p ) const;

  /// Compute distance to given coordinates
  real dist( real const * x ) const;

  /// Compute norm of point representing a vector from the origin
  inline real norm() const;

  /// Return strong relative comparison
  bool operator%( real a );

  /// Compute cross product embedded in R^3 with given vector
  inline point< 3 > cross( point const & p ) const;

  /// Compute dot product with given vector
  inline real dot( point const & p ) const;

  /// Display info
  void disp() const;

private:
  real x_[D];
};

//-----------------------------------------------------------------------------
template < uint D >
inline void point< D >::set( real const * x, uint dim )
{
  dolfin_assert( dim <= D );
  std::copy( x, x + dim, x_ );
}

//-----------------------------------------------------------------------------
template < uint D >
inline real & point< D >::operator[]( uint i )
{
  dolfin_assert( i < D );
  return x_[i];
}

//-----------------------------------------------------------------------------
template < uint D >
inline real const & point< D >::operator[]( uint i ) const
{
  dolfin_assert( i < D );
  return x_[i];
}

//-----------------------------------------------------------------------------
template < uint D >
inline point< D >::operator real *()
{
  return &x_[0];
}

//-----------------------------------------------------------------------------
template < uint D >
inline point< D >::operator real const *() const
{
  return &x_[0];
}

//-----------------------------------------------------------------------------
template < uint D >
inline point< D > & point< D >::operator+=( point< D > const & p )
{
  for ( size_t i = 0; i < D; ++i )
  {
    x_[i] += p.x_[i];
  }
  return *this;
}

//-----------------------------------------------------------------------------
template < uint D >
inline point< D > point< D >::operator+( point< D > const & p ) const
{
  return ( point( *this ) += p );
}

//-----------------------------------------------------------------------------
template < uint D >
inline point< D > & point< D >::operator-=( point< D > const & p )
{
  for ( size_t i = 0; i < D; ++i )
  {
    x_[i] -= p.x_[i];
  }
  return *this;
}

//-----------------------------------------------------------------------------
template < uint D >
inline point< D > point< D >::operator-( point< D > const & p ) const
{
  return ( point( *this ) -= p );
}

//-----------------------------------------------------------------------------
template < uint D >
inline point< D > & point< D >::operator*=( real a )
{
  for ( size_t i = 0; i < D; ++i )
  {
    x_[i] *= a;
  }
  return *this;
}

//-----------------------------------------------------------------------------
template < uint D >
inline point< D > point< D >::operator*( real a ) const
{
  return ( point< D >( *this ) *= a );
}

//-----------------------------------------------------------------------------
template < uint D >
inline point< D > & point< D >::operator/=( real a )
{
  for ( size_t i = 0; i < D; ++i )
  {
    x_[i] /= a;
  }
  return *this;
}

//-----------------------------------------------------------------------------
template < uint D >
inline point< D > point< D >::operator/( real a ) const
{
  return ( point< D >( *this ) /= a );
}

//-----------------------------------------------------------------------------
template < uint D >
inline real point< D >::dist( point< D > const & p ) const
{
  real d = 0.;
  for ( uint i = 0; i < D; ++i )
  {
    d += ( p.x_[i] - x_[i] ) * ( p.x_[i] - x_[i] );
  }
  return std::sqrt( d );
}

//-----------------------------------------------------------------------------
template < uint D >
inline real point< D >::dist( real const * x ) const
{
  real d = 0.;
  for ( uint i = 0; i < D; ++i )
  {
    d += ( x[i] - x_[i] ) * ( x[i] - x_[i] );
  }
  return std::sqrt( d );
}

//-----------------------------------------------------------------------------
template < uint D >
inline real point< D >::norm() const
{
  real d = 0.;
  for ( uint i = 0; i < D; ++i )
  {
    d += x_[i] * x_[i];
  }
  return std::sqrt( d );
}

//-----------------------------------------------------------------------------
template < uint D >
inline bool point< D >::operator%( real a )
{
  return srelcmp( norm(), a, DOLFIN_EPS );
}

//-----------------------------------------------------------------------------
template < uint D >
inline real point< D >::dot( point< D > const & p ) const
{
  real d = 0.;
  for ( uint i = 0; i < D; ++i )
  {
    d += x_[i] * p.x_[i];
  }
  return d;
}

//-----------------------------------------------------------------------------
template < uint D >
inline void point< D >::disp() const
{
  section( "point" );
  message(
    "( %+e, %+e, %+e )", x_[0], ( D > 1 ? x_[1] : 0 ), ( D > 2 ? x_[2] : 0 ) );
  end();
}

//-----------------------------------------------------------------------------
template <>
inline point< 3 > point< 1 >::cross( point const & ) const
{
  return point< 3 >();
}

//-----------------------------------------------------------------------------
template <>
inline point< 3 > point< 2 >::cross( point const & p ) const
{
  return point< 3 >( 0., 0., x_[0] * p.x_[1] - x_[1] * p.x_[0] );
}

//-----------------------------------------------------------------------------
template <>
inline point< 3 > point< 3 >::cross( point const & p ) const
{
  return point< 3 >( x_[1] * p.x_[2] - x_[2] * p.x_[1],
                     x_[2] * p.x_[0] - x_[0] * p.x_[2],
                     x_[0] * p.x_[1] - x_[1] * p.x_[0] );
}

//-----------------------------------------------------------------------------
/// Determinant R^2
inline real operator^( point< 2 > const & p, point< 2 > const & q )
{
  return p[0] * q[1] - p[1] * q[0];
}

//-----------------------------------------------------------------------------
/// Determinant R^3
inline point< 3 > operator^( point< 3 > const & p, point< 3 > const & q )
{
  return p.cross( q );
}

//-----------------------------------------------------------------------------
template < uint D >
inline point< D > operator*( real a, point< D > const & p )
{
  return p * a;
}

//-----------------------------------------------------------------------------
template < uint D >
inline LogStream & operator<<( LogStream & ss, point< D > const & p )
{
  ss << "[ point x = ( " << p[0];
  for ( uint i = 1; i < D; ++i )
    ss << ", " << p[i];
  ss << " ) ]";
  return ss;
}

//-----------------------------------------------------------------------------
template < uint D >
inline std::ostream & operator<<( std::ostream & ss, point< D > const & p )
{
  ss << "[ point x = ( " << p[0];
  for ( uint i = 1; i < D; ++i )
    ss << ", " << p[i];
  ss << " ) ]";
  return ss;
}

//-----------------------------------------------------------------------------

using Point = class point<>;

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_POINT_H */
