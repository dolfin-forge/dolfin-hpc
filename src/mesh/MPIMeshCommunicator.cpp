// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Magnus Vikstrøm, 2007.
// Modified by Anders Logg, 2007.
// Modified by Niclas Jansson, 2008-2015.
// Modified by Balthasar Reuter, 2013.
//
// First added:  2007-05-30
// Last changed: 2015-01-31

#include <dolfin/mesh/MPIMeshCommunicator.h>

#include <dolfin/common/timing.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MPIMeshCommunicator::MPIMeshCommunicator()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MPIMeshCommunicator::~MPIMeshCommunicator()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(Mesh& mesh,
                                     MeshFunction<uint>& distribution)
{
  uint const ddim = distribution.dim();
  if (ddim == 0)
  {
    distributeVertices(mesh, distribution);
  }
  else if (ddim == mesh.topology().dim())
  {
    distributeCells(mesh, distribution);
  }
  else
  {
    error("MPIMeshCommunicator : unimplemented for dimension %u", ddim);
  }
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distributeVertices(Mesh& mesh,
                                             MeshFunction<uint>& distribution)
{
  if (!mesh.is_distributed())
  {
    return;
  }

#if HAVE_MPI

  message(1, "MPIMeshCommunicator : distribute vertices");
  tic();

  uint const rank = MPI::processNumber();
  uint const pe_size = MPI::numProcesses();
  MeshTopology& topology = mesh.topology();
  uint const tdim = topology.dim();
  MeshGeometry& geometry = mesh.geometry();
  uint const gdim = geometry.dim();

  // Save global num vertices to check consistency
  dolfin_assert(topology.distdata()[0].is_finalized());
  uint const num_global_vertices = topology.global_size(0);

  //
  dolfin_assert(topology.entities_exist(0));
  if (topology.size(0) != distribution.size())
  {
    error("MPIMeshCommunicator : mismatch between number of vertices and size "
          "of the distribution");
  }
  for (uint d = 1; d <= tdim; ++d)
  {
    if (topology.entities_exist(d))
    {
      error("MPIMeshCommunicator : distribution by vertices but entities of "
            "dimension %u exist",
            d);
    }
  }

  DistributedData distdata1;
  Array<real> coords;
  Array<real> *sendbuf_x = new Array<real> [pe_size];
  Array<uint> *sendbuf_i = new Array<uint> [pe_size];

  // Process mesh entities according to distribution
  uint index = 0;
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    if (v->is_owned())
    {
      uint const owner = distribution.get(*v);
      if (owner != rank)
      {
        sendbuf_i[owner].push_back(v->global_index());
        for (uint d = 0; d < gdim; ++d)
        {
          sendbuf_x[owner].push_back(v->x()[d]);
        }
      }
      else
      {
        distdata1.set_map(index, v->global_index());
        ++index;
        for (uint d = 0; d < gdim; ++d)
        {
          coords.push_back(v->x()[d]);
        }
      }
    }
  }

  // Exchange the vertices
  MPI_Status status;
  uint src;
  uint dest;
  uint maxsend_i = sendbuf_i[0].size();
  for (uint j = 1; j < pe_size; ++j)
  {
    maxsend_i = std::max(maxsend_i, (uint) sendbuf_i[j].size());
  }
  uint maxrecv_i = 0;
  MPI::numGlobalSum(maxsend_i, maxrecv_i);
  uint * recvbuf_i = new uint[maxrecv_i];
  uint maxrecv_x = maxrecv_i * gdim;
  real * recvbuf_x = new real[maxrecv_x];
  int recv_count;
  for (uint j = 1; j < pe_size; ++j)
  {
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&sendbuf_i[dest][0], sendbuf_i[dest].size(), MPI_UNSIGNED,
                 dest, 1, &recvbuf_i[0], maxrecv_i, MPI_UNSIGNED, src, 1,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
    MPI_Sendrecv(&sendbuf_x[dest][0], sendbuf_x[dest].size(), MPI_DOUBLE, dest,
                 2, &recvbuf_x[0], recv_count, MPI_DOUBLE, src, 2,
                 MPI::DOLFIN_COMM, &status);

    for (int k = 0; k < recv_count; ++k)
    {
      distdata1.set_map(index, recvbuf_i[k]);
      ++index;
      for (uint d = 0; d < gdim; ++d)
      {
        coords.push_back(recvbuf_x[k * gdim + d]);
      }
    }
  }

  // Cleanup and finalize distributed data
  delete[] recvbuf_x;
  delete[] recvbuf_i;
  delete[] sendbuf_x;
  delete[] sendbuf_i;
  distdata1.finalize();

  // Clear mesh then update mesh topology and geometry
  // Do not use the mesh editor since it clears all data including periodic
  // constraints: redistribution of the mesh would destroy all this information
  uint const num_local_vertices = distdata1.local_size();
  topology.init(tdim);
  topology.init(0 , num_local_vertices);
  topology.distdata()[0] = distdata1;
  topology.finalize();
  if(num_global_vertices == topology.global_size(0))
  {
    error("MPIMeshCommunicator : invalid global number of vertices after "
          "distribution by vertices produced an invalid %u != %u",
          num_global_vertices, topology.global_size(0));
  }
  geometry.init(gdim, num_local_vertices);
  geometry.set(coords);
  geometry.finalize();

  //
  tocd(1);

#endif /* HAVE_MPI */

}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distributeCells(Mesh& mesh,
                                          MeshFunction<uint>& distribution)
{
  if (!mesh.is_distributed())
  {
    return;
  }

#if HAVE_MPI

  message(1, "MPIMeshCommunicator : distribute cells");
  tic();

  error("Unimplemented");

  uint const rank = MPI::processNumber();
  uint const pe_size = MPI::numProcesses();
  MeshTopology& topology = mesh.topology();
  uint const tdim = topology.dim();
  MeshGeometry& geometry = mesh.geometry();
  uint const gdim = geometry.dim();


  //
  tocd(1);

#endif /* HAVE_MPI */

}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
