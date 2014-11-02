// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2006.
// Modified by Aurélien Larcher, 2014.
//
// First added:  2006-06-12
// Last changed: 2014-12-01

#ifndef __POINT_H
#define __POINT_H

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/types.h>

namespace dolfin
{

/// A Point represents a point in R^3 with coordinates x, y, z, or,
/// alternatively, a vector in R^3, supporting standard operations
/// like the norm, distances, scalar and vector products etc.

class Point
{
public:

  static uint const max_size = 3;

  /// Create a point at (x, y, z)
  Point(const real x = 0.0, const real y = 0.0, const real z = 0.0);

  /// Copy constructor
  Point(Point const& p);

  /// Destructor
  ~Point();

  /// Return coordinate in direction i
  real& operator[](uint i);

  /// Return coordinate in direction i
  real const& operator[](uint i) const;

  /// Return x-coordinate
  real x() const;

  /// Return y-coordinate
  real y() const;

  /// Return z-coordinate
  real z() const;

  /// Compute sum of two points
  Point operator+(Point const& p) const;

  /// Compute difference of two points
  Point operator-(Point const& p) const;

  /// Add given point
  Point const& operator+=(Point const& p);

  /// Subtract given point
  Point const& operator-=(Point const& p);

  /// Multiplication with scalar
  Point operator*(real a) const;

  /// Incremental multiplication with scalar
  Point const& operator*=(real a);

  /// Division by scalar
  Point operator/(real a) const;

  /// Incremental division by scalar
  Point const& operator/=(real a);

  /// Assignment operator
  Point const& operator=(Point const& p);

  /// Compute distance to given point
  real distance(Point const& p) const;

  /// Compute norm of point representing a vector from the origin
  real norm() const;

  /// Compute cross product with given vector
  const Point cross(Point const& p) const;

  /// Compute dot product with given vector
  real dot(Point const& p) const;

  /// Output
  friend LogStream& operator<<(LogStream& stream, Point const& p);

  /// Display info
  void disp() const;

private:

  real x_[Point::max_size];

};

//--- INLINES -----------------------------------------------------------------

inline real& Point::operator[](uint i)
{
  dolfin_assert(i < Point::max_size);
  return x_[i];
}

//-----------------------------------------------------------------------------
inline real const& Point::operator[](uint i) const
{
  dolfin_assert(i < Point::max_size);
  return x_[i];
}

//-----------------------------------------------------------------------------
inline real Point::x() const
{
  return x_[0];
}

//-----------------------------------------------------------------------------
inline real Point::y() const
{
  return x_[1];
}

//-----------------------------------------------------------------------------
inline real Point::z() const
{
  return x_[2];
}

//-----------------------------------------------------------------------------
inline Point Point::operator+(Point const& p) const
{
  Point q(x_[0] + p.x_[0], x_[1] + p.x_[1], x_[2] + p.x_[2]);
  return q;
}

//-----------------------------------------------------------------------------
inline Point Point::operator-(Point const& p) const
{
  Point q(x_[0] - p.x_[0], x_[1] - p.x_[1], x_[2] - p.x_[2]);
  return q;
}

//-----------------------------------------------------------------------------
inline Point const& Point::operator+=(Point const& p)
{
  x_[0] += p.x_[0];
  x_[1] += p.x_[1];
  x_[2] += p.x_[2];
  return *this;
}

//-----------------------------------------------------------------------------
inline Point const& Point::operator-=(Point const& p)
{
  x_[0] -= p.x_[0];
  x_[1] -= p.x_[1];
  x_[2] -= p.x_[2];
  return *this;
}

//-----------------------------------------------------------------------------
inline Point Point::operator*(real a) const
{
  Point p(a * x_[0], a * x_[1], a * x_[2]);
  return p;
}

//-----------------------------------------------------------------------------
inline Point operator*(real a, Point const& p)
{
  return p * a;
}

//-----------------------------------------------------------------------------
inline Point const& Point::operator*=(real a)
{
  x_[0] *= a;
  x_[1] *= a;
  x_[2] *= a;
  return *this;
}

//-----------------------------------------------------------------------------
inline Point Point::operator/(real a) const
{
  Point p(x_[0] / a, x_[1] / a, x_[2] / a);
  return p;
}

//-----------------------------------------------------------------------------
inline Point const& Point::operator/=(real a)
{
  x_[0] /= a;
  x_[1] /= a;
  x_[2] /= a;
  return *this;
}

//-----------------------------------------------------------------------------
inline Point const& Point::operator=(Point const& p)
{
  x_[0] = p.x_[0];
  x_[1] = p.x_[1];
  x_[2] = p.x_[2];
  return *this;
}

//-----------------------------------------------------------------------------

}

#endif
