// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/mesh/RefinementPattern.h>

#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
uint RefinementPattern::num_refined_vertices(Mesh const& mesh) const
{
  uint val = 0;
  for (uint i = 0; i <= mesh.topology().dim(); ++i)
  {
    val += this->num_refined_vertices(i) * mesh.size(i);
  }
  return val;
}
//-----------------------------------------------------------------------------
uint RefinementPattern::num_refined_cells(Mesh const& mesh) const
{
  return (this->num_refined_cells() * mesh.num_cells());
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
