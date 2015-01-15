// Copyright (C) 2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-11-23
// Last changed: 2007-11-23

#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/UnitInterval.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UnitInterval::UnitInterval(uint nx) :
    Mesh()
{
  if (nx < 1) error("Size of unit interval must be at least 1.");

  rename("mesh", "Mesh of the unit interval (0,1)");

  // Open mesh for editing
  MeshEditor editor(*this, CellType::interval, 1, 1);

  // Create vertices and cells:
  editor.initVertices((nx + 1));
  editor.initCells(nx);

  // Create main vertices:
  for (uint ix = 0; ix <= nx; ++ix)
  {
    real const x = static_cast<real>(ix) / static_cast<real>(nx);
    editor.addVertex(ix, x);
  }

  // Create intervals
  for (uint ix = 0; ix < nx; ++ix)
  {
    uint const connectivity[2] = { ix, ix + 1 };
    editor.addCell(ix, &connectivity[0]);
  }

  // Close mesh editor
  editor.close();

}
//-----------------------------------------------------------------------------

}
