// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2007.
//
// First added:  2005-12-02
// Last changed: 2007-1-06

#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/UnitSquare.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UnitSquare::UnitSquare(uint nx, uint ny, Type type) :
    Mesh()
{

  if (nx < 1 || ny < 1)
  {
    error("Size of unit square must be at least 1 in each dimension.");
  }

  rename("mesh", "Mesh of the unit square (0,1) x (0,1)");

  // Open mesh for editing
  MeshEditor editor(*this, CellType::triangle, 2, 2);

  // Create vertices and cells:
  if (type == crisscross)
  {
    editor.initVertices((nx + 1) * (ny + 1) + nx * ny);
    editor.initCells(4 * nx * ny);
  }
  else
  {
    editor.initVertices((nx + 1) * (ny + 1));
    editor.initCells(2 * nx * ny);
  }

  // Create main vertices:
  uint vertex = 0;
  real x[2] = { 0.0 };
  for (uint iy = 0; iy <= ny; ++iy)
  {
    x[1] = static_cast<real>(iy) / static_cast<real>(ny);
    for (uint ix = 0; ix <= nx; ++ix)
    {
      x[0] = static_cast<real>(ix) / static_cast<real>(nx);
      editor.addVertex(vertex++, x);
    }
  }

  // Create midpoint vertices if the mesh type is crisscross
  if (type == crisscross)
  {
    for (uint iy = 0; iy < ny; iy++)
    {
      x[1] = (static_cast<real>(iy) + 0.5) / static_cast<real>(ny);
      for (uint ix = 0; ix < nx; ix++)
      {
        x[0] = (static_cast<real>(ix) + 0.5) / static_cast<real>(nx);
        editor.addVertex(vertex++, x);
      }
    }
  }

  // Create triangles
  uint cell = 0;
  if (type == crisscross)
  {
    for (uint iy = 0; iy < ny; iy++)
    {
      for (uint ix = 0; ix < nx; ix++)
      {
        uint const v0 = iy * (nx + 1) + ix;
        uint const v1 = v0 + 1;
        uint const v2 = v0 + (nx + 1);
        uint const v3 = v1 + (nx + 1);
        uint const vmid = (nx + 1) * (ny + 1) + iy * nx + ix;

        // Note that v0 < v1 < v2 < v3 < vmid.
        uint const connectivity[12] = { v0, v1, vmid, v0, v2, vmid, v1, v3,
                                        vmid, v2, v3, vmid };

        editor.addCell(cell++, &connectivity[0]);
        editor.addCell(cell++, &connectivity[3]);
        editor.addCell(cell++, &connectivity[6]);
        editor.addCell(cell++, &connectivity[9]);
      }
    }
  }
  else if (type == left)
  {
    for (uint iy = 0; iy < ny; iy++)
    {
      for (uint ix = 0; ix < nx; ix++)
      {
        uint const v0 = iy * (nx + 1) + ix;
        uint const v1 = v0 + 1;
        uint const v2 = v0 + (nx + 1);
        uint const v3 = v1 + (nx + 1);

        uint const connectivity[6] = { v0, v1, v2, v1, v2, v3 };

        editor.addCell(cell++, &connectivity[0]);
        editor.addCell(cell++, &connectivity[3]);
      }
    }
  }
  else
  {
    for (uint iy = 0; iy < ny; iy++)
    {
      for (uint ix = 0; ix < nx; ix++)
      {
        uint const v0 = iy * (nx + 1) + ix;
        uint const v1 = v0 + 1;
        uint const v2 = v0 + (nx + 1);
        uint const v3 = v1 + (nx + 1);

        uint const connectivity[6] = { v0, v1, v3, v0, v2, v3 };

        editor.addCell(cell++, &connectivity[0]);
        editor.addCell(cell++, &connectivity[3]);
      }
    }
  }

  // Close mesh editor
  editor.close();

}
//-----------------------------------------------------------------------------

}
