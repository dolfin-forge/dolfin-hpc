// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
// Modified by Balthasar Reuter, 2013
//

#include <dolfin/mesh/DVertex.h>

using namespace dolfin;
//------------------------------------------------------------------------------
#ifdef HAVE_LIBGEOM
DVertex::DVertex() : id(0), glb_id(-1), patch_id(-1), u(-1.f), v(-1.f),
         cells(0), p(0.0, 0.0, 0.0), deleted(false), on_boundary(false), 
         shared(false), ghosted(false), owner(-1)
#else
DVertex::DVertex() : id(0), glb_id(-1), cells(0), p(0.0, 0.0, 0.0), 
         deleted(false), on_boundary(false), shared(false),
         ghosted(false), owner(-1)				
#endif
{
}
//------------------------------------------------------------------------------
