// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman 2006.
// Modified by Aurelien Larcher 2015.
//
// First added:  2006-06-02
// Last changed: 2015-02-09

#include <dolfin/mesh/Edge.h>

#include <cmath>

namespace dolfin
{

//-----------------------------------------------------------------------------
real Edge::length() const
{
  uint const * vertices = entities(0);
  dolfin_assert(vertices);
  MeshGeometry const& geom = mesh_.geometry();
  real const * p0 = geom.x(vertices[0]);
  real const * p1 = geom.x(vertices[1]);
  real l = 0;
  for (uint n = 0; n < gdim_; ++n)
  {
    l += (p1[n] - p0[n]) * (p1[n] - p0[n]);
  }
  return std::sqrt(l);
}
//-----------------------------------------------------------------------------
Point Edge::midpoint() const
{
  uint const * vertices = entities(0);
  dolfin_assert(vertices);
  MeshGeometry const& geom = mesh_.geometry();
  real const * p0 = geom.x(vertices[0]);
  real const * p1 = geom.x(vertices[1]);
  Point p;
  for (uint n = 0; n < gdim_; ++n)
  {
    p[n] = 0.5 * (p0[n] + p1[n]);
  }
  return p;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
