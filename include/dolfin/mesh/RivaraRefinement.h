// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_RIVARA_REFINEMENT_H
#define __DOLFIN_RIVARA_REFINEMENT_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

class Mesh;

namespace RivaraRefinement
{

/// Refine simplicial mesh locally by recursive edge bisection
void refine( Mesh& mesh, MeshValues< bool, Cell > & cell_marker,
             real tf = 0.0, real tb = 0.0, real ts = 0.0,
             bool balance = true );

} // end namespce RivaraRefinement

} // end namespce dolfin

#endif
