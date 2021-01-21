// Copyright (C) 2008 Nuno Lopes.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UNIT_SPHERE_H
#define __DOLFIN_UNIT_SPHERE_H

#include "Mesh.h"

namespace dolfin
{

/// Triangular mesh of the 3D unit sphere.
///
/// Given the number of cells (nx, ny, nz) in each direction,
/// the total number of tetrahedra will be 6*nx*ny*nz and the
/// total number of vertices will be (nx + 1)*(ny + 1)*(nz + 1).

class UnitSphere : public Mesh
{
public:

  UnitSphere(uint nx);

private:

  auto transformx(real x, real y, real z) -> real;
  auto transformy(real x, real y, real z) -> real;
  auto transformz(real x, real y, real z) -> real;
  auto max(real x, real y, real z) -> real;

};

}

#endif
