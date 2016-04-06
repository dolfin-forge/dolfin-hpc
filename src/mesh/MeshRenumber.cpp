// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/mesh/MeshRenumber.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
bool MeshRenumber::renumber(Mesh& mesh)
{
  if (!mesh.is_distributed())
  {
    return false;
  }

  bool renumbered = false;

#ifdef HAVE_MPI

  MeshTopology& topology = mesh.topology();
  MeshDistributedData& distdata = topology.distdata();

  uint const tdim = topology.dim();
  uint const rank = MPI::processNumber();
  uint const pe_size = MPI::numProcesses();

  //---------------------------------------------------------------------------

  /*
   * Renumber vertices: compacting global numbering per rank.
   *
   */

  if (topology.entities_exist(0) && !distdata[0].valid_numbering)
  {
    DistributedData& distdata0 = distdata[0];
    dolfin_assert(distdata0.is_finalized());
    distdata0.renumber_global();
    //
    distdata0.valid_numbering = true;
  }

  //---------------------------------------------------------------------------

  /*
   * Renumber edges/faces.
   *
   */

  for (uint d = 0; d < tdim; ++d)
  {
    if (topology.entities_exist(d) && !distdata[d].valid_numbering)
    {
      DistributedData& distdata0 = distdata[d];
      //
      distdata0.valid_numbering = true;
    }
  }

  //---------------------------------------------------------------------------

  /*
   * Renumber cells: compacting global numbering per rank.
   *
   */

  if (topology.entities_exist(tdim) && !distdata[tdim].valid_numbering)
  {
    DistributedData& distdata0 = distdata[tdim];
    dolfin_assert(distdata0.is_finalized());
    distdata0.renumber_global();
    //
    distdata0.valid_numbering = true;
  }

#endif /* HAVE_MPI */

  return renumbered;
}

} /* namespace dolfin */

