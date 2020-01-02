// Copyright (C) 2008 Solveig Bruvoll and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __MESH_ALE_H
#define __MESH_ALE_H

#include <dolfin/common/types.h>

namespace dolfin
{

class Mesh;
class BoundaryMesh;

/// This file provides functionality useful for implementation of
/// ALE (Arbitrary Lagrangian-Eulerian) methods, in particular
/// moving the boundary vertices of a mesh and then interpolating
/// the new coordinates for the interior vertices accordingly.

namespace ALE
{

/// List of available methods for ALE mesh movement
enum ALEType
{
	lagrange,
	hermite,
	harmonic,
	elastic
};

/// Move coordinates of mesh according to new boundary coordinates
void move( Mesh & mesh, BoundaryMesh & new_boundary, ALEType type = lagrange );

} // namespace ALE

} // namespace dolfin

#endif // __MESH_ALE_H
