// Copyright (C) 2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-11-23
// Last changed: 2007-11-23

#include <dolfin/mesh/UnitInterval.h>

#include <dolfin/math/LinearDistribution.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/MeshEditor.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UnitInterval::UnitInterval(uint nx) :
    Mesh()
{
  if (nx < 1) error("Size of unit interval must be at least 1.");

  rename("mesh", "Mesh of the unit interval (0,1)");

  // Open mesh for editing
  MeshEditor editor(*this, CellType::interval, 1);

  //
  uint const rank = MPI::processNumber();
  uint const pe_size = MPI::numProcesses();
  LinearDistribution cdist(nx, pe_size);

  // Create vertices and cells:
  editor.init_vertices(cdist.size + 1, cdist.global_size + 1);
  editor.init_cells(cdist.size, cdist.global_size);

  // Create main vertices:
  for (uint ix = 0; ix <= cdist.size; ++ix)
  {
    real const x = real(cdist.offset + ix) / real(cdist.global_size);
    editor.add_vertex(ix, &x);
  }

  // Create intervals
  for (uint ix = 0; ix < cdist.size; ++ix)
  {
    uint const connectivity[2] = { cdist.offset + ix, cdist.offset + ix + 1 };
    editor.add_cell(ix, &connectivity[0]);
  }

  if (is_distributed())
  {
    for (uint ix = 0; ix <= cdist.size; ++ix)
    {
      distdata()[0].set_map(ix, cdist.offset + ix);
    }
    if (rank > 0)
    {
      distdata()[0].set_ghost(cdist.offset, rank - 1);
    }
    if (rank < pe_size - 1)
    {
      distdata()[0].set_shared_adj(cdist.offset + cdist.size, rank + 1);
    }
  }

  // Close mesh editor
  editor.close();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
