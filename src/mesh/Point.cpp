// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2006.
// Modified by Aurélien Larcher, 2014.
//
// First added:  2006-06-12
// Last changed: 2014-06-12

#include <dolfin/mesh/Point.h>

#include <cmath>

namespace dolfin
{

//-----------------------------------------------------------------------------
Point::Point(const real x, const real y, const real z)
{
  _x[0] = x;
  _x[1] = y;
  _x[2] = z;
}
//-----------------------------------------------------------------------------
Point::Point(Point const& p)
{
  _x[0] = p._x[0];
  _x[1] = p._x[1];
  _x[2] = p._x[2];
}
//-----------------------------------------------------------------------------
Point::~Point()
{
}
//-----------------------------------------------------------------------------
real Point::distance(Point const& p) const
{
  const real dx = p._x[0] - _x[0];
  const real dy = p._x[1] - _x[1];
  const real dz = p._x[2] - _x[2];

  return std::sqrt(dx * dx + dy * dy + dz * dz);
}
//-----------------------------------------------------------------------------
real Point::norm() const
{
  return std::sqrt(_x[0] * _x[0] + _x[1] * _x[1] + _x[2] * _x[2]);
}
//-----------------------------------------------------------------------------
const Point Point::cross(Point const& p) const
{
  Point q;

  q._x[0] = _x[1] * p._x[2] - _x[2] * p._x[1];
  q._x[1] = _x[2] * p._x[0] - _x[0] * p._x[2];
  q._x[2] = _x[0] * p._x[1] - _x[1] * p._x[0];

  return q;
}
//-----------------------------------------------------------------------------
real Point::dot(Point const& p) const
{
  return _x[0] * p._x[0] + _x[1] * p._x[1] + _x[2] * p._x[2];
}
//-----------------------------------------------------------------------------
LogStream& operator<<(LogStream& stream, Point const& p)
{
  stream << "[ Point x = " << p.x() << " y = " << p.y() << " z = " << p.z()
         << " ]";
  return stream;
}
//-----------------------------------------------------------------------------

}
