// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007.
//
// First added:  2005-12-02
// Last changed: 2007-12-06

#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/UnitCube.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UnitCube::UnitCube(uint nx, uint ny, uint nz) :
    Mesh()
{

  if (nx < 1 || ny < 1 || nz < 1)
  {
    error("Size of unit cube must be at least 1 in each dimension.");
  }

  rename("mesh", "Mesh of the unit cube (0,1) x (0,1) x (0,1)");

  // Open mesh for editing
  MeshEditor editor(*this, CellType::tetrahedron, 3, 3);

  // Create vertices
  editor.initVertices((nx + 1) * (ny + 1) * (nz + 1));
  uint vertex = 0;
  real x[3] = { 0.0 };
  for (uint iz = 0; iz <= nz; iz++)
  {
    x[2] = static_cast<real>(iz) / static_cast<real>(nz);
    for (uint iy = 0; iy <= ny; iy++)
    {
      x[1] = static_cast<real>(iy) / static_cast<real>(ny);
      for (uint ix = 0; ix <= nx; ix++)
      {
        x[0] = static_cast<real>(ix) / static_cast<real>(nx);
        editor.addVertex(vertex++, x);
      }
    }
  }

  // Create tetrahedra
  editor.initCells(6 * nx * ny * nz);
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

        editor.addCell(cell++, &connectivity[0]);
        editor.addCell(cell++, &connectivity[4]);
        editor.addCell(cell++, &connectivity[8]);
        editor.addCell(cell++, &connectivity[12]);
        editor.addCell(cell++, &connectivity[16]);
        editor.addCell(cell++, &connectivity[20]);
      }
    }
  }

  // Close mesh editor
  editor.close();

}
//-----------------------------------------------------------------------------

}
