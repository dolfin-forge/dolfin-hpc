// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
// Modified by Balthasar Reuter, 2013
//

#include <dolfin/mesh/DVertex.h>

using namespace dolfin;
//------------------------------------------------------------------------------
DVertex::DVertex() : id(0), glb_id(-1), cells(0), p(0.0, 0.0, 0.0), 
         deleted(false), on_boundary(false), shared(false),
         ghosted(false), owner(-1)
{
}
//------------------------------------------------------------------------------