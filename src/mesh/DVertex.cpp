// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
// Modified by Balthasar Reuter, 2013
//

#include <dolfin/mesh/DVertex.h>

namespace dolfin
{

//------------------------------------------------------------------------------
DVertex::DVertex() :
    id(0),
    glb_id(-1),
    cells(0),
    p(0.0, 0.0, 0.0),
    deleted(false),
    shared(false),
    ghosted(false),
    owner(-1)
{
}
//------------------------------------------------------------------------------
DVertex::DVertex(Vertex const& v) :
    id(v.index()),
    glb_id(v.global_index()),
    cells(0),
    p(v.point()),
    deleted(false),
    shared(v.is_shared()),
    ghosted(v.is_ghost()),
    owner(v.owner())
{
  if (this->shared)
  {
    this->shared_adj = v.mesh().distdata()[0].get_shared_adj(id);
  }
}
//------------------------------------------------------------------------------

} /* namespace dolfin */
