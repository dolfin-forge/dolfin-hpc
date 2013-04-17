// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
// Modified by Balthasar Reuter, 2013
//

#include <dolfin/mesh/DCell.h>
#include <dolfin/mesh/DVertex.h>

using namespace dolfin;

//------------------------------------------------------------------------------
DCell::DCell() : id(0), parent_id(0), vertices(0), deleted(false), nref(0)
{
}
//-----------------------------------------------------------------------------
bool DCell::has_edge(DVertex *v1, DVertex *v2)
{
  uint found = 0;
  for ( std::vector<DVertex*>::iterator it = vertices.begin() ;
        it != vertices.end(); ++it )
    if ( *it == v1 || *it == v2 )
      found++;
  return (found == 2);
}
//------------------------------------------------------------------------------