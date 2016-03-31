// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007.
// Modified by Nuno Lopes, 2008
//
// First added:  2005-12-02
// Last changed: 2008-06-20

#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/UnitSphere.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UnitSphere::UnitSphere(uint nx) :
    Mesh()
{

  message("UnitSphere is Experimental: It could have a bad quality mesh");

  uint ny = nx;
  uint nz = nx;

  if (nx < 1 || ny < 1 || nz < 1)
  {
    error("Size of unit cube must be at least 1 in each dimension.");
  }

  rename("mesh", "Mesh of the unit cube (0,1) x (0,1) x (0,1)");

  // Open mesh for editing
  MeshEditor editor(*this, CellType::tetrahedron, 3);

  // Create vertices
  editor.init_vertices((nx + 1) * (ny + 1) * (nz + 1));
  uint vertex = 0;
  real trns_x[3] = { 0.0 };
  for (uint iz = 0; iz <= nz; iz++)
  {
    real const z = -1.0 + static_cast<real>(iz) * 2.0 / static_cast<real>(nz);
    for (uint iy = 0; iy <= ny; iy++)
    {
      real const y = -1.0 + static_cast<real>(iy) * 2.0 / static_cast<real>(ny);
      for (uint ix = 0; ix <= nx; ix++)
      {
        real const x = -1.0
            + static_cast<real>(ix) * 2.0 / static_cast<real>(nx);
        trns_x[0] = transformx(x, y, z);
        trns_x[1] = transformy(x, y, z);
        trns_x[2] = transformz(x, y, z);
        editor.add_vertex(vertex++, trns_x);
      }
    }
  }

  // Create tetrahedra
  editor.init_cells(6 * nx * ny * nz);
  uint cell = 0;
  for (uint iz = 0; iz < nz; iz++)
  {
    for (uint iy = 0; iy < ny; iy++)
    {
      for (uint ix = 0; ix < nx; ix++)
      {
        uint const v0 = iz * (nx + 1) * (ny + 1) + iy * (nx + 1) + ix;
        uint const v1 = v0 + 1;
        uint const v2 = v0 + (nx + 1);
        uint const v3 = v1 + (nx + 1);
        uint const v4 = v0 + (nx + 1) * (ny + 1);
        uint const v5 = v1 + (nx + 1) * (ny + 1);
        uint const v6 = v2 + (nx + 1) * (ny + 1);
        uint const v7 = v3 + (nx + 1) * (ny + 1);

        uint const connectivity[24] = { v0, v1, v3, v7, v0, v1, v7, v5, v0, v5,
                                        v7, v4, v0, v3, v2, v7, v0, v6, v4, v7,
                                        v0, v2, v6, v7 };

        editor.add_cell(cell++, &connectivity[0]);
        editor.add_cell(cell++, &connectivity[4]);
        editor.add_cell(cell++, &connectivity[8]);
        editor.add_cell(cell++, &connectivity[12]);
        editor.add_cell(cell++, &connectivity[16]);
        editor.add_cell(cell++, &connectivity[20]);
      }
    }
  }

  // Close mesh editor
  editor.close();

}
//-----------------------------------------------------------------------------
real UnitSphere::transformx(real x, real y, real z)
{
  real retrn = 0.0;
  if (x || y || z)

  {
    retrn = x * max(fabs(x), fabs(y), fabs(z)) / sqrt(x * x + y * y + z * z);
  }
  else
  {
    retrn = x;
  }
  return retrn;
}
//-----------------------------------------------------------------------------
real UnitSphere::transformy(real x, real y, real z)
{
  real retrn = 0.0;
  if (x || y || z)
  {
    retrn = y * max(fabs(x), fabs(y), fabs(z)) / sqrt(x * x + y * y + z * z);
  }
  else
  {
    retrn = y;
  }
  return retrn;
}
//-----------------------------------------------------------------------------
real UnitSphere::transformz(real x, real y, real z)
{
  real retrn = 0.0;
  //maxn transformation
  if (x || y || z)
  {
    retrn = z * max(fabs(x), fabs(y), fabs(z)) / sqrt(x * x + y * y + z * z);
  }
  else
  {
    retrn = z;
  }
  return retrn;
}
//-----------------------------------------------------------------------------
real UnitSphere::max(real x, real y, real z)
{
  real rtrn = 0.0;

  if ((x >= y) * (x >= z))
  {
    rtrn = x;
  }
  else if ((y >= x) * (y >= z))
  {
    rtrn = y;
  }
  else
  {
    rtrn = z;
  }
  return rtrn;
}
//-----------------------------------------------------------------------------

}
