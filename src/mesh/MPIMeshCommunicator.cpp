// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Magnus Vikstrøm, 2007.
// Modified by Anders Logg, 2007.
// Modified by Niclas Jansson, 2008-2015.
// Modified by Balthasar Reuter, 2013.
// Modified by Aurelien Larcher, 2015-2016.
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

#include <algorithm>

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
void MPIMeshCommunicator::distribute(Mesh& mesh, MeshFunction<uint>& dist)
{
  uint const ddim = dist.dim();
  if (ddim == 0)
  {
    distributeVertices(mesh, dist);
  }
  else if (ddim == mesh.topology().dim())
  {
    distributeCells(mesh, dist);
  }
  else
  {
    error("MPIMeshCommunicator : unimplemented for dimension %u", ddim);
  }
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distributeVertices(Mesh& mesh, MeshFunction<uint>& dist)
{
  if (!mesh.is_distributed())
  {
    return;
  }

#if HAVE_MPI

  message(1, "MPIMeshCommunicator : distribute vertices");
  tic();

  uint const rank = MPI::rank();
  uint const pe_size = MPI::size();
  MeshTopology& topology = mesh.topology();
  uint const tdim = topology.dim();
  MeshGeometry& geometry = mesh.geometry();
  uint const gdim = geometry.dim();

  // Save global number of vertices to check consistency
  dolfin_assert(topology.entities_exist(0));
  dolfin_assert(topology.distdata()[0].is_finalized());
  uint const num_global_vertices = topology.global_size(0);
  if (topology.size(0) != dist.size())
  {
    error("MPIMeshCommunicator : mismatch between number of vertices and size "
          "of the distribution");
  }
  for (uint d = 1; d <= tdim; ++d)
  {
    if (topology.entities_exist(d) && (topology.size(d) > 0))
    {
      error("MPIMeshCommunicator : distribution by vertices but entities of "
            "dimension %u exist", d);
    }
  }

  DistributedData distdata1;
  Array<real> coords;
  Array<uint> * sendbuf_v = new Array<uint> [pe_size];
  Array<real> * sendbuf_x = new Array<real> [pe_size];

  // Collect mesh entities according to distribution
  uint vindex = 0;
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    if (v->is_owned())
    {
      uint const owner = dist.get(*v);
      if (owner == rank)
      {
        distdata1.set_map(vindex, v->global_index());
        ++vindex;
        for (uint d = 0; d < gdim; ++d)
        {
          coords.push_back(v->x()[d]);
        }
      }
      else
      {
        sendbuf_v[owner].push_back(v->global_index());
        for (uint d = 0; d < gdim; ++d)
        {
          sendbuf_x[owner].push_back(v->x()[d]);
        }
      }
    }
  }

  // Clear topology and geometry
  topology.clear();
  geometry.clear();

  // Exchange the vertices
  MPI_Status status;
  uint src;
  uint dst;
  uint send_size;
  uint recvmax_v;
  uint recvmax_x;
  for (uint j = 0; j < pe_size; ++j)
  {
    send_size = sendbuf_v[j].size();
    MPI_Reduce(&send_size, &recvmax_v, 1, MPI_UNSIGNED, MPI_MAX, j,
               MPI::DOLFIN_COMM);
    send_size = sendbuf_x[j].size();
    MPI_Reduce(&send_size, &recvmax_x, 1, MPI_UNSIGNED, MPI_SUM, j,
               MPI::DOLFIN_COMM);
  }
  dolfin_assert(recvmax_v > 0);
  // Allocate vertex indices buffer
  uint * recvbuf_v = new uint[recvmax_v];
  // Resize vertex coordinates array to fit new cells
  uint const coords_size = coords.size();
  coords.resize(coords_size + recvmax_x);
  real * recvbuf_x = &coords[coords_size];
  int recv_count;
  for (uint j = 1; j < pe_size; ++j)
  {
    src = (rank - j + pe_size) % pe_size;
    dst = (rank + j) % pe_size;

    // Vertices
    MPI_Sendrecv(&sendbuf_v[dst][0], sendbuf_v[dst].size(), MPI_UNSIGNED,
                 dst, 0, &recvbuf_v[0], recvmax_v, MPI_UNSIGNED, src, 0,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int k = 0; k < recv_count; ++k)
    {
      uint const global_index = recvbuf_v[k];
      if(distdata1.has_global(global_index))
      {
        error("MPIMeshCommunicator : receiving global vertex %u twice",
              global_index);
      }
      distdata1.set_map(vindex, global_index);
      ++vindex;
    }

    // Coordinates
    MPI_Sendrecv(&sendbuf_x[dst][0], sendbuf_x[dst].size(), MPI_DOUBLE, dst, 1,
                 recvbuf_x, recvmax_x, MPI_DOUBLE, src, 1, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
    recvbuf_x += recv_count;
    recvmax_x -= recv_count;

  }

  // Cleanup and finalize distributed data
  delete[] recvbuf_v;
  delete[] sendbuf_x;
  delete[] sendbuf_v;

  // Finalize distributed data
  distdata1.finalize();

  // Update topology
  dolfin_assert(vindex == distdata1.local_size());
  topology.init(tdim);
  topology.init(0 , vindex);
  topology.distdata()[0] = distdata1;
  topology.finalize();
  if(num_global_vertices != topology.global_size(0))
  {
    error("MPIMeshCommunicator : invalid global number of vertices %u != %u",
          num_global_vertices, topology.global_size(0));
  }

  // Update geometry
  dolfin_assert(vindex * gdim == coords.size());
  geometry.init(gdim, vindex);
  geometry.set(coords);
  geometry.finalize();

  //
  tocd(1);

#endif /* HAVE_MPI */

}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distributeCells(Mesh& mesh, MeshFunction<uint>& dist)
{
  if (!mesh.is_distributed())
  {
    return;
  }

#if HAVE_MPI

  message(1, "MPIMeshCommunicator : distribute cells");
  tic();

  uint const rank = MPI::rank();
  uint const pe_size = MPI::size();
  MeshTopology& topology = mesh.topology();
  uint const tdim = topology.dim();
  MeshGeometry& geometry = mesh.geometry();
  uint const gdim = geometry.dim();

  // Save global number of vertices and cells to check consistency
  dolfin_assert(topology.entities_exist(0));
  dolfin_assert(topology.distdata()[0].is_finalized());
  uint const num_global_vertices = topology.global_size(0);
  dolfin_assert(topology.entities_exist(tdim));
  dolfin_assert(topology.distdata()[tdim].is_finalized());
  uint const num_global_cells = topology.global_size(tdim);
  if (topology.size(tdim) != dist.size())
  {
    error("MPIMeshCommunicator : mismatch between number of cells and size of "
          "the distribution");
  }

  DistributedData distdata1;
  Array<real> coords;
  Array<uint> * sendbuf_c = new Array<uint> [pe_size];
  Array<uint> * sendbuf_v = new Array<uint> [pe_size];
  Array<real> * sendbuf_x = new Array<real> [pe_size];

  // Collect mesh entities according to distribution
  Array<uint> cells;
  uint vindex = 0;
  bool * vertex_used = new bool[topology.size(0)];
  std::fill_n(vertex_used, topology.size(0), false);
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    uint const owner = dist.get(*c);
    if (owner == rank)
    {
      for (VertexIterator v(*c); !v.end(); ++v)
      {
        cells.push_back(v->global_index());
        if (!vertex_used[v->index()] && v->is_owned())
        {
          vertex_used[v->index()] = true;
          distdata1.set_map(vindex, v->global_index());
          ++vindex;
          for (uint d = 0; d < gdim; ++d)
          {
            coords.push_back(v->x()[d]);
          }
        }
      }
    }
    else
    {
      for (VertexIterator v(*c); !v.end(); ++v)
      {
        sendbuf_c[owner].push_back(v->global_index());
        if (!vertex_used[v->index()] && v->is_owned())
        {
          vertex_used[v->index()] = true;
          sendbuf_v[owner].push_back(v->global_index());
          for (uint d = 0; d < gdim; ++d)
          {
            sendbuf_x[owner].push_back(v->x()[d]);
          }
        }
      }
    }
  }
  delete [] vertex_used;

  // Clear topology and geometry
  topology.clear();
  geometry.clear();

  // Exchange the processed entities
  MPI_Status status;
  uint src;
  uint dst;
  uint send_size;
  uint recvmax_c;
  uint recvmax_v;
  uint recvmax_x;
  for (uint j = 0; j < pe_size; ++j)
  {
    send_size = sendbuf_c[j].size();
    MPI_Reduce(&send_size, &recvmax_c, 1, MPI_UNSIGNED, MPI_SUM, j,
               MPI::DOLFIN_COMM);
    send_size = sendbuf_v[j].size();
    MPI_Reduce(&send_size, &recvmax_v, 1, MPI_UNSIGNED, MPI_MAX, j,
               MPI::DOLFIN_COMM);
    send_size = sendbuf_x[j].size();
    MPI_Reduce(&send_size, &recvmax_x, 1, MPI_UNSIGNED, MPI_SUM, j,
               MPI::DOLFIN_COMM);
  }
  // Resize cell vertices array to fit new cells
  uint const cells_size = cells.size();
  cells.resize(cells_size + recvmax_c);
  uint * recvbuf_c = &cells[cells_size];
  // Allocate vertex indices buffer
  dolfin_assert(recvmax_v > 0);
  uint * recvbuf_v = new uint[recvmax_v];
  // Resize vertex coordinates array to fit new cells
  uint const coords_size = coords.size();
  coords.resize(coords_size + recvmax_x);
  real * recvbuf_x = &coords[coords_size];
  int recv_count;
  for (uint j = 1; j < pe_size; ++j)
  {
    src = (rank - j + pe_size) % pe_size;
    dst = (rank + j) % pe_size;

    // Cells
    MPI_Sendrecv(&sendbuf_c[dst][0], sendbuf_c[dst].size(), MPI_UNSIGNED, dst,
                 0, recvbuf_c, recvmax_c, MPI_UNSIGNED, src, 0,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
    recvbuf_c += recv_count;
    recvmax_c -= recv_count;

    // Vertices
    MPI_Sendrecv(&sendbuf_v[dst][0], sendbuf_v[dst].size(), MPI_UNSIGNED, dst,
                 1, recvbuf_v, recvmax_v, MPI_UNSIGNED, src, 1,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int k = 0; k < recv_count; ++k)
    {
      uint const global_index = recvbuf_v[k];
      if(distdata1.has_global(global_index))
      {
        error("MPIMeshCommunicator : receiving global vertex %u twice",
              global_index);
      }
      distdata1.set_map(vindex, global_index);
      ++vindex;
    }

    // Coordinates
    MPI_Sendrecv(&sendbuf_x[dst][0], sendbuf_x[dst].size(), MPI_DOUBLE, dst, 2,
                 recvbuf_x, recvmax_x, MPI_DOUBLE, src, 2, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
    recvbuf_x += recv_count;
    recvmax_x -= recv_count;

  }

  // Cleanup buffers
  delete[] recvbuf_v;
  delete[] sendbuf_x;
  delete[] sendbuf_v;
  delete[] sendbuf_c;

  // NOTE: This implementation only works for homogeneous topologies
  //       Check cell data size just in case.
  if ((cells.size() % mesh.type().num_entities(0)) > 0)
  {
    error("MPIMeshCommunicator : inconsistent size of cell buffer '%u'",
          cells.size());
  }
  uint cindex = cells.size() / mesh.type().num_entities(0);

  // Loop over cells as a list of global vertices and determine if vertices are
  // local or not: overwrite the array with local indices to avoid copy.
  // Fill buffer with ghost vertices.
  Array<uint> sendbuf_gv;
  _set<uint> global_gv;
  for (Array<uint>::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    uint const global_index = (*it);
    if (distdata1.has_global(global_index))
    {
      (*it) = distdata1.get_local(global_index);
    }
    else
    {
      (*it) = vindex;
      // Map new vertex
      distdata1.set_map(vindex, global_index);
      ++vindex;
      // Padding for coordinates
      for (uint d = 0; d < gdim; ++d)
      {
        coords.push_back(0.0);
      }
      // Add vertex to ghost buffer
      sendbuf_gv.push_back(global_index);
      global_gv.insert(global_index);
    }
  }

  // Exchange ghost vertices
  uint sendcnt_gv = sendbuf_gv.size();
  uint sendmax_gv = 0;
  MPI::allReduceSum(sendcnt_gv, sendmax_gv);
  dolfin_assert(sendmax_gv > 0);
  uint * sendbck_gv = new uint[sendmax_gv];
  real * sendbck_gx = new real[sendmax_gv * gdim];
  uint * recvbuf_gv = new uint[sendcnt_gv];
  real * recvbuf_gx = new real[sendcnt_gv * gdim];
  for (uint j = 1; j < pe_size; ++j)
  {
    src = (rank - j + pe_size) % pe_size;
    dst = (rank + j) % pe_size;

    // Send ghost vertices to request coordinates
    MPI_Sendrecv(&sendbuf_gv[0], sendbuf_gv.size(), MPI_UNSIGNED, dst, 0,
                 sendbck_gv, sendmax_gv, MPI_UNSIGNED, src, 0, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    uint send_count = 0;
    for (int k = 0; k <recv_count; ++k)
    {
      uint const global_index = sendbck_gv[k];
      if (distdata1.has_global(global_index))
      {
        uint const local_index = distdata1.get_local(global_index);
        // Set source rank as shared adjacent
        distdata1.set_shared_adj(local_index, src);
        // If vertex is owned shared then send coordinates back
        if (!global_gv.count(global_index))
        {
          sendbck_gv[send_count] = global_index;
          std::copy(&coords[local_index * gdim],
                    &coords[local_index * gdim] + gdim,
                    &sendbck_gx[send_count * gdim]);
          ++send_count;
        }
      }
    }

    // Send coordinates back
    MPI_Sendrecv(&sendbck_gv[0], send_count, MPI_UNSIGNED, src, 1,
                 &recvbuf_gv[0], sendcnt_gv, MPI_UNSIGNED, dst, 1,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    sendcnt_gv -= recv_count;

    MPI_Sendrecv(&sendbck_gx[0], send_count * gdim, MPI_DOUBLE, src, 2,
                 &recvbuf_gx[0], recv_count * gdim, MPI_DOUBLE, dst, 2,
                 MPI::DOLFIN_COMM, &status);

    for (int k = 0; k < recv_count; ++k)
    {
      uint const local_index = distdata1.get_local(recvbuf_gv[k]);
      // Set vertex owner and copy coordinates
      distdata1.set_ghost(local_index, dst);
      std::copy(&recvbuf_gx[k * gdim], &recvbuf_gx[k * gdim] + gdim,
                &coords[local_index * gdim]);
    }

  }

  // Cleanup
  global_gv.clear();
  delete [] recvbuf_gx;
  delete [] recvbuf_gv;
  delete [] sendbck_gx;
  delete [] sendbck_gv;

  // Finalize distributed data
  distdata1.finalize();

  // Update topology
  dolfin_assert(vindex == distdata1.local_size());
  topology.init(tdim);
  topology.init(0 , vindex);
  topology.distdata()[0] = distdata1;
  topology.init(tdim , cindex);
  message("Init cells");
  topology(tdim, 0).set(cells);
  topology.finalize();
  if(num_global_vertices != topology.global_size(0))
  {
    error("MPIMeshCommunicator : invalid global number of vertices %u != %u",
          num_global_vertices, topology.global_size(0));
  }
  if(num_global_cells != topology.global_size(tdim))
  {
    error("MPIMeshCommunicator : invalid global number of cells %u != %u",
          num_global_cells, topology.global_size(tdim));
  }

  // Update geometry
  dolfin_assert(vindex * gdim == coords.size());
  geometry.init(gdim, vindex);
  geometry.set(coords);
  geometry.finalize();

  //
  tocd(1);

#endif /* HAVE_MPI */

}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
