// Copyright (C) 2013 Balthasar Reuter.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_BARYCENTRIC_POINT_H
#define __DOLFIN_BARYCENTRIC_POINT_H

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace dolfin
{

/// A barycentric point represents a point with respect to a barycentric
/// coordinate system.

class BarycentricPoint
{
public:
  /// Create a barycentric point with coordinates
  /// (lambda_x, lambda_y, lambda_z)
  BarycentricPoint( real lambda_x = 0.,
                    real lambda_y = 0.,
                    real lambda_z = 0. );

  /// Copy constructor
  BarycentricPoint( BarycentricPoint const & p );

  /// Destructor
  ~BarycentricPoint()
  {
  }

  /// Return reference to coordinate in direction i
  inline real & operator[]( uint i );

  /// Return coordinate in direction i
  inline real operator[]( uint i ) const;

  /// Compute transformation matrix (as row major) and translation vector
  /// with respect to a cell
  static void getTransformation( Cell const & c,
                                 real ** matrix, real ** vector );

  /// Compute euclidian point to barycentric coordinate with respect to a
  /// certain cell. Transformation vector and matrix can be supplied to avoid
  /// recomputation.
  Point point( Cell const & c,
               real * matrix = NULL, real * vector = NULL ) const;

private:
  real _x[3];
};

//-----------------------------------------------------------------------------
inline real & BarycentricPoint::operator[]( uint i )
{
  dolfin_assert( i < 3 );
  return _x[i];
}

//-----------------------------------------------------------------------------
inline real BarycentricPoint::operator[]( uint i ) const
{
  dolfin_assert( i < 3 );
  return _x[i];
}

//-----------------------------------------------------------------------------
inline void BarycentricPoint::getTransformation( Cell const & c,
                                                 real ** matrix,
                                                 real ** vector )
{
  uint         n     = c.num_entities( 0 ) - 1;
  uint         dim   = c.mesh().geometry_dimension();
  const uint * verts = c.entities( 0 );
  real *       mat   = new real[dim * n];
  real *       vec   = new real[dim];
  std::copy( c.mesh().geometry().x( verts[n] ),
             c.mesh().geometry().x( verts[n] ) + dim,
             vec );

  for ( uint i( 0 ); i < n; ++i )
  {
    Point column = c.mesh().geometry().point( verts[i] )
                   - c.mesh().geometry().point( verts[n] );
    for ( uint row( 0 ); row < dim; ++row )
      mat[row * n + i] = column[row];
  }

  *matrix = mat;
  *vector = vec;
}

//-----------------------------------------------------------------------------
inline Point BarycentricPoint::point( Cell const & c,
                                      real * matrix,
                                      real * vector ) const
{
  if ( matrix == 0 || vector == 0 )
    getTransformation( c, &matrix, &vector );

  uint  n   = c.num_entities( 0 ) - 1;
  uint  dim = c.mesh().geometry_dimension();
  Point p;

  for ( uint row( 0 ); row < dim; ++row )
  {
    p[row] = vector[row];
    for ( uint col( 0 ); col < n; ++col )
      p[row] += matrix[row * n + col] * _x[col];
  }

  return p;
}

//-----------------------------------------------------------------------------

} // namespace

#endif
