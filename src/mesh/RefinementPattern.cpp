// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/RefinementPattern.h>

#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

auto RefinementPattern::num_refined_vertices( Mesh const & mesh ) const
  -> size_t
{
  size_t val = 0;
  for ( size_t i = 0; i <= mesh.topology_dimension(); ++i )
  {
    val += this->num_refined_vertices( i ) * mesh.size( i );
  }
  return val;
}

//-----------------------------------------------------------------------------

auto RefinementPattern::num_refined_cells( Mesh const & mesh ) const -> size_t
{
  return ( this->num_refined_cells() * mesh.num_cells() );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
