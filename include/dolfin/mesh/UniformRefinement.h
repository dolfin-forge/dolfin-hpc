// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_UNIFORM_REFINEMENT_H
#define __DOLFIN_MESH_UNIFORM_REFINEMENT_H

namespace dolfin
{

class Mesh;

/// This file implements uniform mesh refinement for different mesh types.

namespace UniformRefinement
{

/// Refine mesh uniformly according to default pattern
void refine( Mesh& mesh );

} // namespace UniformRefinement

} // namespace dolfin

#endif /* __DOLFIN_MESH_UNIFORM_REFINEMENT_H */
