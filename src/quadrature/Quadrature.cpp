// Copyright (C) 2003-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/log/dolfin_log.h>
#include <dolfin/quadrature/Quadrature.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
Quadrature::Quadrature(unsigned int n)
{
  this->n = n;

  points = new real[n];
  weights = new real[n];

  for (unsigned int i = 0; i < n; i++)
    weights[i] = 0;

  m = 1.0;
}
//-----------------------------------------------------------------------------
Quadrature::~Quadrature()
{
  if ( points )
    delete [] points;
  points = nullptr;

  if ( weights )
    delete [] weights;
  weights = nullptr;
}
//-----------------------------------------------------------------------------
auto Quadrature::size() const -> unsigned int
{
  return n;
}
//-----------------------------------------------------------------------------
auto Quadrature::point(unsigned int i) const -> real
{
  dolfin_assert(i < n);
  return points[i];
}
//-----------------------------------------------------------------------------
auto Quadrature::weight(unsigned int i) const -> real
{
  dolfin_assert(i < n);
  return weights[i];
}
//-----------------------------------------------------------------------------
auto Quadrature::measure() const -> real
{
  return m;
}
//-----------------------------------------------------------------------------
