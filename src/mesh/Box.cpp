// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007.
// Modified by Nuno Lopes, 2008.
//
// First added:  2005-12-02
// Last changed: 2008-06-19

#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Box.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Box::Box(real a, real b, real c, real d, real e, real f, uint nx, uint ny,
         uint nz) :
    Mesh()
{

  if (nx < 1 || ny < 1 || nz < 1)
  {
    error("Size of box must be at least 1 in each dimension.");
  }

  rename("mesh", "Mesh of the cuboid (a,b) x (c,d) x (e,f)");

  // Open mesh for editing
  MeshEditor editor(*this, CellType::tetrahedron, 3, 3);

  // Create vertices
  editor.initVertices((nx + 1) * (ny + 1) * (nz + 1));
  uint vertex = 0;
  real x[3] = { 0.0 };
  for (uint iz = 0; iz <= nz; ++iz)
  {
    x[2] = e + (static_cast<real>(iz)) * (f - e) / static_cast<real>(nz);
    for (uint iy = 0; iy <= ny; ++iy)
    {
      x[1] = c + (static_cast<real>(iy)) * (d - c) / static_cast<real>(ny);
      for (uint ix = 0; ix <= nx; ++ix)
      {
        x[0] = a + (static_cast<real>(ix)) * (b - a) / static_cast<real>(nx);
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
