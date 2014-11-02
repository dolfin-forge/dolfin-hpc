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
  x_[0] = x;
  x_[1] = y;
  x_[2] = z;
}
//-----------------------------------------------------------------------------
Point::Point(Point const& p)
{
  x_[0] = p.x_[0];
  x_[1] = p.x_[1];
  x_[2] = p.x_[2];
}
//-----------------------------------------------------------------------------
Point::~Point()
{
}
//-----------------------------------------------------------------------------
real Point::distance(Point const& p) const
{
  const real dx = p.x_[0] - x_[0];
  const real dy = p.x_[1] - x_[1];
  const real dz = p.x_[2] - x_[2];

  return std::sqrt(dx * dx + dy * dy + dz * dz);
}
//-----------------------------------------------------------------------------
real Point::norm() const
{
  return std::sqrt(x_[0] * x_[0] + x_[1] * x_[1] + x_[2] * x_[2]);
}
//-----------------------------------------------------------------------------
const Point Point::cross(Point const& p) const
{
  Point q;

  q.x_[0] = x_[1] * p.x_[2] - x_[2] * p.x_[1];
  q.x_[1] = x_[2] * p.x_[0] - x_[0] * p.x_[2];
  q.x_[2] = x_[0] * p.x_[1] - x_[1] * p.x_[0];

  return q;
}
//-----------------------------------------------------------------------------
real Point::dot(Point const& p) const
{
  return x_[0] * p.x_[0] + x_[1] * p.x_[1] + x_[2] * p.x_[2];
}
//-----------------------------------------------------------------------------
LogStream& operator<<(LogStream& stream, Point const& p)
{
  stream << "[ Point x = " << p.x() << " y = " << p.y() << " z = " << p.z()
         << " ]";
  return stream;
}
//-----------------------------------------------------------------------------
void Point::disp() const
{
  // Begin indentation
  cout << "Point" << endl;
  begin("-----");
  cout << "( " << x_[0] << ", " << x_[1] << ", " << x_[2] << " )" << endl;
  // End indentation
  end();
}
//-----------------------------------------------------------------------------

}
