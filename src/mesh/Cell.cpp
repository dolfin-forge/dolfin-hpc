// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman 2006.
//
// First added:  2006-01-01
// Last changed: 2007-04-16

#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Point Cell::midpoint() const
{
  MeshGeometry const& geometry = this->mesh().geometry();
  uint const* vertices = this->entities(0);
  uint const num_vertices = this->numEntities(0);
  Point p;
  for (uint v = 0; v < num_vertices; ++v)
  {
    real const* x = geometry.x(vertices[v]);
    for (uint i = 0; i < geometry.dim(); ++i)
    {
      p[i] += x[i];
    }
  }
  p /= real(num_vertices);
  return p;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
