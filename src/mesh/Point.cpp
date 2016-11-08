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
Point::Point(real const * x, uint gdim)
{
  std::fill(x_, x_ + Point::MAX_SIZE, 0.0);
  std::copy(x, x + gdim, x_);
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
  return std::sqrt((p.x_[0] - x_[0]) * (p.x_[0] - x_[0]) +
                   (p.x_[1] - x_[1]) * (p.x_[1] - x_[1]) +
                   (p.x_[2] - x_[2]) * (p.x_[2] - x_[2]));
}
//-----------------------------------------------------------------------------
real Point::norm() const
{
  return std::sqrt(x_[0] * x_[0] + x_[1] * x_[1] + x_[2] * x_[2]);
}
//-----------------------------------------------------------------------------
Point Point::cross(Point const& p) const
{
  return Point(x_[1] * p.x_[2] - x_[2] * p.x_[1],
               x_[2] * p.x_[0] - x_[0] * p.x_[2],
               x_[0] * p.x_[1] - x_[1] * p.x_[0]);
}
//-----------------------------------------------------------------------------
real Point::dot(Point const& p) const
{
  return x_[0] * p.x_[0] + x_[1] * p.x_[1] + x_[2] * p.x_[2];
}
//-----------------------------------------------------------------------------
LogStream& operator<<(LogStream& ss, Point const& p)
{
  ss << "[ Point x = " << p.x() << " y = " << p.y() << " z = " << p.z() << " ]";
  return ss;
}
//-----------------------------------------------------------------------------
void Point::disp() const
{
  section("Point");
  message("( %+e, %+e, %+e )", x_[0], x_[1], x_[2]);
  endblock();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
