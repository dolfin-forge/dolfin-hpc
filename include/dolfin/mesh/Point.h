// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_POINT_H
#define __DOLFIN_POINT_H

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Space.h>

#include <iostream>

namespace dolfin
{

//-----------------------------------------------------------------------------

template < uint dim >
class point;

using Point = class point< Space::MAX_DIMENSION >;

//-----------------------------------------------------------------------------

/// A point represents a point in R^d with coordinates x, y, z, or,
/// alternatively, a vector in R^d, supporting standard operations
/// like the norm, distances, scalar and vector products etc.

template < uint dim >
class point : public std::array< real, dim >
{
public:
  static uint const MAX_SIZE = dim;

public:
  /// default constructor
  inline point();

  /// Copy constructor
  inline point( point const & other );

  /// Create a point<1> (x)
  explicit inline point( real const x );

  /// Create a point<2> (x, y)
  explicit inline point( real const x, real const y );

  /// Create a point<3> (x, y, z)
  explicit inline point( real const x, real const y, real const z );

  /// Create a point from a given coordinate
  explicit inline point( real const * x );

  /// Destructor
  ~point() = default;

  /// Assignment operator
  inline point & operator=( point const & p );

  /// arithmetic operators
  inline point operator+( point const & p ) const;
  inline point operator-( point const & p ) const;
  inline point operator*( real a ) const;
  inline point operator/( real a ) const;

  inline point & operator+=( point const & p );
  inline point & operator-=( point const & p );
  inline point & operator*=( real a );
  inline point & operator/=( real a );

  inline point & set( real * x );

  /// Compute distance to given point
  inline real dist( point const & p ) const;

  /// Compute distance to given coordinates
  inline real dist( real const * x ) const;

  /// Compute norm of point representing a vector from the origin
  inline real norm() const;

  /// Compute cross product embedded in R^3 with given vector
  inline point cross( point const & p ) const;

  /// Compute cross product embedded in R^3 with given coordinate vector
  inline point cross( real const * x ) const;

  /// Compute dot product with given vector
  inline real dot( point const & p ) const;

  /// Compute dot product with given coordinates
  inline real dot( real const * x ) const;

  /// Display info
  void disp() const;

private:
  using array_type   = typename std::array< real, dim >;
  real * const data_ = this->data();
};

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim >::point()
{
  this->fill( 0.0 );
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim >::point( point const & other )
  : array_type( other )
{
}

//-----------------------------------------------------------------------------

template <>
inline point< 1 >::point( real const x )
{
  data_[0] = x;
}

//-----------------------------------------------------------------------------

template <>
inline point< 2 >::point( real const x, real const y )
{
  data_[0] = x;
  data_[1] = y;
}

//-----------------------------------------------------------------------------

template <>
inline point< 3 >::point( real const x, real const y, real const z )
{
  data_[0] = x;
  data_[1] = y;
  data_[2] = z;
}

//-----------------------------------------------------------------------------

template <>
inline point< 1 >::point( real const * x )
{
  data_[0] = x[0];
}

//-----------------------------------------------------------------------------

template <>
inline point< 2 >::point( real const * x )
{
  data_[0] = x[0];
  data_[1] = x[1];
}

//-----------------------------------------------------------------------------

template <>
inline point< 3 >::point( real const * x )
{
  data_[0] = x[0];
  data_[1] = x[1];
  data_[2] = x[2];
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim > & point< dim >::operator=( point< dim > const & p )
{
  for ( uint i = 0; i < dim; ++i )
    data_[i] = p[i];

  return *this;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim > point< dim >::operator+( point< dim > const & p ) const
{
  point tmp( *this );
  tmp += p;
  return tmp;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim > point< dim >::operator-( point< dim > const & p ) const
{
  point tmp( *this );
  tmp -= p;
  return tmp;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim > point< dim >::operator*( real a ) const
{
  point tmp( *this );
  tmp *= a;
  return tmp;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim > point< dim >::operator/( real a ) const
{
  point tmp( *this );
  tmp /= a;
  return tmp;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim > & point< dim >::operator+=( point< dim > const & p )
{
  for ( uint i = 0; i < dim; ++i )
    data_[i] += p.data_[i];

  return *this;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim > & point< dim >::operator-=( point< dim > const & p )
{
  for ( uint i = 0; i < dim; ++i )
    data_[i] -= p.data_[i];

  return *this;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim > & point< dim >::operator*=( real a )
{
  for ( uint i = 0; i < dim; ++i )
    data_[i] *= a;

  return *this;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim > & point< dim >::operator/=( real a )
{
  for ( uint i = 0; i < dim; ++i )
    data_[i] /= a;

  return *this;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim > & point< dim >::set( real * x )
{
  for ( uint i = 0; i < dim; ++i )
    data_[i] = x[i];

  return *this;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline real point< dim >::dist( point< dim > const & p ) const
{
  real d = 0.;

  for ( uint i = 0; i < dim; ++i )
    d += ( p.data_[i] - data_[i] ) * ( p.data_[i] - data_[i] );

  return std::sqrt( d );
}

//-----------------------------------------------------------------------------

template < uint dim >
inline real point< dim >::dist( real const * x ) const
{
  real d = 0.;

  for ( uint i = 0; i < dim; ++i )
    d += ( x[i] - data_[i] ) * ( x[i] - data_[i] );

  return std::sqrt( d );
}

//-----------------------------------------------------------------------------

template < uint dim >
inline real point< dim >::norm() const
{
  real d = 0.;

  for ( uint i = 0; i < dim; ++i )
    d += data_[i] * data_[i];

  return std::sqrt( d );
}

//-----------------------------------------------------------------------------

template <>
inline point< 3 > point< 3 >::cross( point const & p ) const
{
  return point< 3 >( data_[1] * p.data_[2] - data_[2] * p.data_[1],
                     data_[2] * p.data_[0] - data_[0] * p.data_[2],
                     data_[0] * p.data_[1] - data_[1] * p.data_[0] );
}

//-----------------------------------------------------------------------------

template <>
inline point< 3 > point< 3 >::cross( real const * x ) const
{
  return point< 3 >( data_[1] * x[2] - data_[2] * x[1],
                     data_[2] * x[0] - data_[0] * x[2],
                     data_[0] * x[1] - data_[1] * x[0] );
}

//-----------------------------------------------------------------------------

template < uint dim >
inline real point< dim >::dot( point< dim > const & p ) const
{
  real d = 0.;

  for ( uint i = 0; i < dim; ++i )
    d += data_[i] * p.data_[i];

  return d;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline real point< dim >::dot( real const * x ) const
{
  real d = 0.;

  for ( uint i = 0; i < dim; ++i )
    d += data_[i] * x[i];

  return d;
}

//-----------------------------------------------------------------------------

template < uint D >
inline void point< D >::disp() const
{
  section( "point" );
  message( "( %+e, %+e, %+e )",
           data_[0], ( D > 1 ? data_[1] : 0. ), ( D > 2 ? data_[2] : 0. ) );
  end();
}

//-----------------------------------------------------------------------------

template < uint dim >
inline point< dim > operator*( real a, point< dim > const & p )
{
  return p * a;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline LogStream & operator<<( LogStream & ss, point< dim > const & p )
{
  ss << "[ point x = ( " << p[0];
  for ( uint i = 1; i < dim; ++i )
    ss << ", " << p[i];
  ss << " ) ]";
  return ss;
}

//-----------------------------------------------------------------------------

template < uint dim >
inline std::ostream & operator<<( std::ostream & ss, point< dim > const & p )
{
  ss << "[ point x = ( " << p[0];
  for ( uint i = 1; i < dim; ++i )
    ss << ", " << p[i];
  ss << " ) ]";
  return ss;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_POINT_H */
