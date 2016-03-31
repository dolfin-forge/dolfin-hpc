// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-06-02
// Last changed: 2006-06-02

#include <dolfin/mesh/Facet.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Point Facet::midpoint() const
{
  MeshGeometry const& geometry = this->mesh().geometry();
  uint const* vertices = this->entities(0);
  uint const num_vertices = this->num_entities(0);
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

