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
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(MeshValues<uint, Vertex>& dist)
{

  if (!dist.mesh().is_distributed())
  {
    return;
  }

#if HAVE_MPI

  message(1, "MPIMeshCommunicator : distribute vertices");
  tic();

  Mesh& mesh = dist.mesh();
  uint const pe_rank = MPI::rank();
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

  Array<uint> * sendbuf_v = new Array<uint> [pe_size];
  Array<real> * sendbuf_x = new Array<real> [pe_size];

  // Collect mesh entities according to distribution
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    if (v->is_owned())
    {
      uint const owner = dist(*v);
      sendbuf_x[owner].append(v->x(), v->x() + gdim);
      sendbuf_v[owner].push_back(v->global_index());
    }
  }

  // Swap local indices
  Array<uint> iverts; iverts.swap(sendbuf_v[pe_rank]);
  dolfin_assert(sendbuf_v[pe_rank].size() == 0);

  // Swap local coordinates
  Array<real> coords; coords.swap(sendbuf_x[pe_rank]);
  dolfin_assert(sendbuf_x[pe_rank].size() == 0);

  // Create distributed data
  DistributedData distdata(mesh.distdata()[0].comm());

  // Clear mesh using swap with new instance
  {
    Mesh new_mesh(mesh.type(), mesh.space());
    new_mesh.swap(mesh);
  }
  dolfin_assert(mesh.topology().size(0) == 0);

  // Exchange the vertices
  MPI_Status status;
  uint src;
  uint dst;
  uint recvmax_v;
  for (uint j = 0; j < pe_size; ++j)
  {
    uint s = sendbuf_v[j].size();
    MPI_Reduce(&s, &recvmax_v, 1, MPI_UNSIGNED, MPI_SUM, j, distdata.comm());
  }
  uint recvmax_x = recvmax_v * gdim;
  // Resize vertex indices
  uint const iverts_size = iverts.size();
  iverts.resize(iverts_size + recvmax_v);
  uint * recvbuf_v = iverts.ptr() + iverts_size;
  // Resize vertex coordinates array to fit new cells
  uint const coords_size = coords.size();
  coords.resize(coords_size + recvmax_x);
  real * recvbuf_x = coords.ptr() + coords_size;
  int recv_count;
  for (uint j = 1; j < pe_size; ++j)
  {
    src = (pe_rank - j + pe_size) % pe_size;
    dst = (pe_rank + j) % pe_size;

    // Vertices
    MPI_Sendrecv(&sendbuf_v[dst][0], sendbuf_v[dst].size(), MPI_UNSIGNED,
                 dst, 0, &recvbuf_v[0], recvmax_v, MPI_UNSIGNED, src, 0,
                 distdata.comm(), &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
    recvbuf_v += recv_count;
    recvmax_v -= recv_count;

    // Coordinates
    MPI_Sendrecv(&sendbuf_x[dst][0], sendbuf_x[dst].size(), MPI_DOUBLE, dst, 1,
                 recvbuf_x, recvmax_x, MPI_DOUBLE, src, 1, distdata.comm(),
                 &status);
    MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
    recvbuf_x += recv_count;
    recvmax_x -= recv_count;

  }

  // Cleanup and finalize distributed data
  delete[] sendbuf_x;
  delete[] sendbuf_v;

  // Map local vertex indices
  uint vindex = 0;
  distdata.set_size(iverts.size());
  for (Array<uint>::const_iterator v = iverts.begin(); v != iverts.end(); ++v)
  {
#if DEBUG
    if (distdata.has_global(*v))
    {
      error("MPIMeshCommunicator : duplicate global vertex index %u", *v);
    }
#endif
    distdata.set_map(vindex++, *v);
  }
  dolfin_assert(vindex == distdata.local_size());

  // Finalize distributed data
  distdata.finalize();

  // Update topology
  dolfin_assert(vindex == distdata.local_size());
  topology.set_distributed();
  topology.init(0 , vindex);
  topology.distdata()[0].swap(distdata);
  topology.finalize();
  dolfin_assert(vindex == topology.distdata()[0].local_size());
  if(num_global_vertices != topology.global_size(0))
  {
    error("MPIMeshCommunicator : vertex distribution :\n"
          "invalid global number of vertices %u != %u",
          num_global_vertices, topology.global_size(0));
  }

  // Update geometry
  dolfin_assert(vindex * gdim == coords.size());
  geometry.assign(coords);
  geometry.finalize();

  //
  tocd(1);

#endif /* HAVE_MPI */

}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(MeshValues<uint, Cell>& dist, MeshData * D)
{
  if (!dist.mesh().is_distributed())
  {
    return;
  }

#if HAVE_MPI

  message(1, "MPIMeshCommunicator : distribute cells");
  tic();

  Mesh& mesh = dist.mesh();
  uint const pe_rank = MPI::rank();
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

  Array<uint> * sendbuf_c = new Array<uint> [pe_size];
  Array<uint> * sendbuf_v = new Array<uint> [pe_size];
  Array<real> * sendbuf_x = new Array<real> [pe_size];

  /// Support the same cell and vertex function as before
  uint const numRV = ( D ? D->size<real, Vertex>() : 0 );
  Array<real> * RV = ( numRV ? new Array<real> [pe_size] : NULL );
  uint const numUC = ( D ? D->size<bool, Cell>() + D->size<uint, Cell>() : 0 );
  Array<uint> * UC = ( numUC ? new Array<uint> [pe_size] : NULL );

  /// Check that
  if (D && (D->size() != (numRV + numUC)))
  {
    error("MPIMeshCommunicator : transferring <real, Vertex>, <bool, Cell>, "
          "and <uint, Cell> only is supported.");
  }

  // Collect mesh entities according to distribution
  bool * vertex_used = (topology.size(0) ? new bool[topology.size(0)]() : NULL);
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    uint const owner = dist(*c);
    for (VertexIterator v(*c); !v.end(); ++v)
    {
      sendbuf_c[owner].push_back(v->global_index());
      if (!vertex_used[v->index()])
      {
        vertex_used[v->index()] = true;
        if (v->is_owned())
        {
          sendbuf_v[owner].push_back(v->global_index());
          sendbuf_x[owner].append(v->x(), v->x() + gdim);

          // Transfer vertex functions
          if (RV)
          {
            for (MeshData::iterator<real, Vertex> it(*D); it.valid(); ++it)
            {
              RV[owner].push_back((*it)(*v));
            }
          }
        }
      }
    }

    // Transfer cell functions
    if (UC)
    {
      for (MeshData::iterator<bool, Cell> it(*D); it.valid(); ++it)
      {
        UC[owner].push_back((*it)(*c));
      }
      for (MeshData::iterator<uint, Cell> it(*D); it.valid(); ++it)
      {
        UC[owner].push_back((*it)(*c));
      }
    }
  }
  delete [] vertex_used;

  // Swap local vertex indices
  Array<uint> iverts; iverts.swap(sendbuf_v[pe_rank]);
  dolfin_assert(sendbuf_v[pe_rank].size() == 0);

  // Swap local coordinates
  Array<real> coords; coords.swap(sendbuf_x[pe_rank]);
  dolfin_assert(sendbuf_x[pe_rank].size() == 0);

  // Swap local cells
  Array<uint> vcells; vcells.swap(sendbuf_c[pe_rank]);
  dolfin_assert(sendbuf_c[pe_rank].size() == 0);

  // Create distributed data
  DistributedData distdata(mesh.distdata()[0].comm());

  // Clear mesh using swap with new instance
  {
    Mesh new_mesh(mesh.type(), mesh.space());
    new_mesh.swap(mesh);
  }

  // Exchange the processed entities
  MPI_Status status;
  uint src;
  uint dst;
  uint recvmax[2] = { 0, 0 };
  for (uint j = 0; j < pe_size; ++j)
  {
    uint sendcnt[2] = { (uint)sendbuf_c[j].size(), (uint)sendbuf_v[j].size() };
    MPI_Reduce(sendcnt, recvmax, 2, MPI_UNSIGNED, MPI_SUM, j, distdata.comm());
  }
  uint recvmax_x = recvmax[1] * gdim;
  // Resize cell vertices array to fit new cells
  uint const vcells_size = vcells.size();
  vcells.resize(vcells_size + recvmax[0]);
  uint * recvbuf_c = vcells.ptr() + vcells_size;
  // Resize vertex indices
  uint const iverts_size = iverts.size();
  iverts.resize(iverts_size + recvmax[1]);
  uint * recvbuf_v = iverts.ptr() + iverts_size;
  // Resize vertex coordinates
  uint const coords_size = coords.size();
  coords.resize(coords_size + recvmax_x);
  real * recvbuf_x = coords.ptr() + coords_size;

  // Naive MeshValues exchange until I fix the ghost bug in the template class
  uint recvmaxUC  = recvmax[0] / mesh.type().num_entities(0) * numUC;
  uint *recvbufUC = NULL;
  if (UC)
  {
    uint const sizeUC = UC[pe_rank].size();
    UC[pe_rank].resize(sizeUC + recvmaxUC);
    recvbufUC = &UC[pe_rank][sizeUC];
  }
  uint recvmaxRV  = recvmax[1] * numRV;
  real *recvbufRV = NULL;
  if (RV)
  {
    uint const sizeRV = RV[pe_rank].size();
    RV[pe_rank].resize(sizeRV + recvmaxRV);
    recvbufRV = &RV[pe_rank][sizeRV];
  }

  int recv_count;
  for (uint j = 1; j < pe_size; ++j)
  {
    src = (pe_rank - j + pe_size) % pe_size;
    dst = (pe_rank + j) % pe_size;

    // Cells
    MPI_Sendrecv(&sendbuf_c[dst][0], sendbuf_c[dst].size(), MPI_UNSIGNED, dst,
                 0, recvbuf_c, recvmax[0], MPI_UNSIGNED, src, 0,
                 distdata.comm(), &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
    recvbuf_c += recv_count;
    recvmax[0]-= recv_count;

    // Vertices
    MPI_Sendrecv(&sendbuf_v[dst][0], sendbuf_v[dst].size(), MPI_UNSIGNED, dst,
                 1, recvbuf_v, recvmax[1], MPI_UNSIGNED, src, 1,
                 distdata.comm(), &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
    recvbuf_v += recv_count;
    recvmax[1]-= recv_count;

    // Coordinates
    MPI_Sendrecv(&sendbuf_x[dst][0], sendbuf_x[dst].size(), MPI_DOUBLE, dst, 2,
                 recvbuf_x, recvmax_x, MPI_DOUBLE, src, 2, distdata.comm(),
                 &status);
    MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
    recvbuf_x += recv_count;
    recvmax_x -= recv_count;

    // Transfer cell functions
    if (UC)
    {
      MPI_Sendrecv(&UC[dst][0], UC[dst].size(), MPI_UNSIGNED, dst, 3,
                   recvbufUC  , recvmaxUC     , MPI_UNSIGNED, src, 3,
                   distdata.comm(), &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
      dolfin_assert(recvmaxUC >= (uint) recv_count);
      recvbufUC += recv_count;
      recvmaxUC -= recv_count;
    }

    // Transfer vertex functions
    if (RV)
    {
      MPI_Sendrecv(&RV[dst][0], RV[dst].size(), MPI_DOUBLE, dst, 4,
                   recvbufRV  , recvmaxRV     , MPI_DOUBLE, src, 4,
                   distdata.comm(), &status);
      MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
      dolfin_assert(recvmaxRV >= (uint) recv_count);
      recvbufRV += recv_count;
      recvmaxRV -= recv_count;
    }

  }

  // Swap local mesh values to new Array.
  Array<uint> mUC; if (UC) { mUC.swap(UC[pe_rank]); delete [] UC; }
  Array<real> mRV; if (RV) { mRV.swap(RV[pe_rank]); delete [] RV; }

  // Cleanup buffers
  delete[] sendbuf_x;
  delete[] sendbuf_v;
  delete[] sendbuf_c;

  // Map local vertex indices
  uint vindex = 0;
  for (Array<uint>::const_iterator it = iverts.begin(); it != iverts.end(); ++it)
  {
#if DEBUG
    if (distdata.has_global(*it))
    {
      error("MPIMeshCommunicator : duplicate global vertex index %u", *it);
    }
#endif
    distdata.set_map(vindex++, *it);
  }

  // NOTE: This implementation only works for homogeneous topologies
  //       Check cell data size just in case.
  if ((vcells.size() % mesh.type().num_entities(0)) > 0)
  {
    error("MPIMeshCommunicator : inconsistent size of cell buffer '%u'",
          vcells.size());
  }
  uint cindex = vcells.size() / mesh.type().num_entities(0);

  // Loop over cells as a list of global vertices and determine if vertices are
  // local or not: overwrite the array with local indices to avoid copy.
  // Fill buffer with ghost vertices.
  Array<uint> sendbuf_gv;
  _set<uint> global_gv;
  for (Array<uint>::iterator it = vcells.begin(); it != vcells.end(); ++it)
  {
    uint const global_index = (*it);
    if (distdata.has_global(global_index))
    {
      (*it) = distdata.get_local(global_index);
    }
    else
    {
      (*it) = vindex;
      // Map new vertex
      distdata.set_map(vindex, global_index);
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
  MPI::all_reduce<MPI::sum>(sendcnt_gv, sendmax_gv);
  dolfin_assert(sendmax_gv > 0);
  uint * sendbck_gv = new uint[sendmax_gv];
  real * sendbck_gx = new real[sendmax_gv * gdim];
  dolfin_assert(sendcnt_gv > 0);
  uint * recvbuf_gv = new uint[sendcnt_gv];
  real * recvbuf_gx = new real[sendcnt_gv * gdim];
  for (uint j = 1; j < pe_size; ++j)
  {
    src = (pe_rank - j + pe_size) % pe_size;
    dst = (pe_rank + j) % pe_size;

    // Send ghost vertices to request coordinates
    MPI_Sendrecv(&sendbuf_gv[0], sendbuf_gv.size(), MPI_UNSIGNED, dst, 0,
                 sendbck_gv, sendmax_gv, MPI_UNSIGNED, src, 0, distdata.comm(),
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    uint send_count = 0;
    for (int k = 0; k <recv_count; ++k)
    {
      uint const global_index = sendbck_gv[k];
      if (distdata.has_global(global_index))
      {
        uint const local_index = distdata.get_local(global_index);
        // Set source rank as shared adjacent
        distdata.set_shared_adj(local_index, src);
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
                 distdata.comm(), &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    sendcnt_gv -= recv_count;

    MPI_Sendrecv(&sendbck_gx[0], send_count * gdim, MPI_DOUBLE, src, 2,
                 &recvbuf_gx[0], recv_count * gdim, MPI_DOUBLE, dst, 2,
                 distdata.comm(), &status);

    for (int k = 0; k < recv_count; ++k)
    {
      uint const local_index = distdata.get_local(recvbuf_gv[k]);
      // Set vertex owner and copy coordinates
      distdata.set_ghost(local_index, dst);
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
  distdata.finalize();

  // Update topology
  dolfin_assert(vindex == distdata.local_size());
  topology.set_distributed();
  topology.init(0 , vindex);
  topology.distdata()[0].swap(distdata);
  topology.init(tdim , cindex);
  topology(tdim, 0).set(vcells);
  topology.finalize();
  dolfin_assert(vindex == topology.distdata()[0].local_size());
  if(num_global_vertices != topology.global_size(0))
  {
    error("MPIMeshCommunicator : cell distribution :\n"
          "invalid global number of vertices %u != %u",
          num_global_vertices, topology.global_size(0));
  }
  if(num_global_cells != topology.global_size(tdim))
  {
    error("MPIMeshCommunicator : cell distribution :\n"
          "invalid global number of cells %u != %u",
          num_global_cells, topology.global_size(tdim));
  }

  // Update geometry
  dolfin_assert(vindex * gdim == coords.size());
  geometry.assign(coords);
  geometry.finalize();

  // Recreate mesh functions
  if (UC)
  {
    uint ii = 0;
    for (MeshData::iterator<bool, Cell> it(*D); it.valid(); ++it, ++ii)
    {
      MeshValues<bool, Cell> M(mesh);
      uint const nUC = mUC.size();
      dolfin_assert(nUC == M.size() * numUC);
      for (uint j = ii, k = 0; j < nUC; j+=numUC, ++k) { M(k) = mUC[j]; }
      (*it).swap(M);
    }
    for (MeshData::iterator<uint, Cell> it(*D); it.valid(); ++it, ++ii)
    {
      MeshValues<uint, Cell> M(mesh);
      uint const nUC = mUC.size();
      dolfin_assert(nUC == M.size() * numUC);
      for (uint j = ii, k = 0; j < nUC; j+=numUC, ++k) { M(k) = mUC[j]; }
      (*it).swap(M);
    }
  }
  if (RV)
  {
    uint ii = 0;
    for (MeshData::iterator<real, Vertex> it(*D); it.valid(); ++it, ++ii)
    {
      MeshValues<real, Vertex> M(mesh);
      uint const nRV = mRV.size();
      dolfin_assert(nRV == (M.size()- mesh.distdata()[0].num_ghost()) * numRV);
      for (uint j = ii, k = 0; j < nRV; j+=numRV, ++k) { M(k) = mRV[j]; }
      (*it).swap(M);
    }
  }

  //
  tocd(1);

#endif /* HAVE_MPI */

}
//-----------------------------------------------------------------------------
template<class E>
void MPIMeshCommunicator::check(Mesh& mesh)
{
  if (!mesh.is_distributed()) return;

  uint const tdim = mesh.topology().dim();
  uint const edim = entity_dimension<E>(mesh);

  if (edim > tdim)
  {
    error("MPIMeshCommunicator : invalid entity dimension %u", edim);
  }

#if HAVE_MPI

  uint const pe_rank = dolfin::MPI::rank();
  uint const pe_size = dolfin::MPI::size();
  DistributedData& dist = mesh.distdata()[edim];

  message("MPIMeshCommunicator : check distribution for dimension %u", edim);

  // Check shared entities adjacency
  {
    Array<uint> * sbuf = new Array<uint> [pe_size];
    uint e_count = 0;
    for (typename E::shared e(mesh); !e.end(); ++e, ++e_count)
    {
      e.adj_enqueue(sbuf, e.global_index());
      // Check that entity adjacency is a subset of adjacent ranks
      for (_set<uint>::const_iterator it = e.adj().begin(); it != e.adj().end(); ++it)
      {
        if (dist.get_adj_ranks().count(*it) == 0)
        {
          error("Invalid adjacent rank for entity %u", e.index());
        }
        dolfin_assert(sbuf[*it].back() == e.global_index());
      }
      // Check that owned entities indices are within process range
      if (e.is_owned() && !dist.in_range(e.global_index()))
      {
        error("Global index of owned entity %u is not in range", e.index());
      }
    }
    // Check that shared entities indices were counted correctly
    dolfin_assert(e_count == mesh.topology().num_shared(edim));

    // Swap local entities to array
    Array<uint> rbuf; rbuf.swap(sbuf[pe_rank]);
    dolfin_assert(sbuf[pe_rank].size() == 0);

    // Exchange entities
    MPI_Status status;
    uint src;
    uint dst;
    uint recv_max;
    for (uint j = 0; j < pe_size; ++j)
    {
      uint s = sbuf[j].size();
      MPI_Reduce(&s, &recv_max, 1, MPI_UNSIGNED, MPI_SUM, j, dist.comm());
      // Check that no duplicate global entity was added
      std::set<uint> global_indices(sbuf[j].begin(), sbuf[j].end());
      if (global_indices.size() != sbuf[j].size())
      {
        error("Duplicate global indices for entities of dimension %u", edim);
      }
    }
    uint const rbuf_size = rbuf.size();
    rbuf.resize(rbuf_size + recv_max);
    uint * recv_buf = rbuf.ptr() + rbuf_size;
    int recv_count;
    for (uint j = 1; j < pe_size; ++j)
    {
      src = (pe_rank - j + pe_size) % pe_size;
      dst = (pe_rank + j) % pe_size;

      MPI_Sendrecv(&sbuf[dst][0], sbuf[dst].size(), MPI_UNSIGNED, dst, 0,
                   &recv_buf[0] , recv_max        , MPI_UNSIGNED, src, 0,
                   dist.comm(), &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

      for (int k = 0; k < recv_count; ++k)
      {
        uint const gindex = recv_buf[k];
        // Check that receive entity global index exists on current rank
        if (dist.has_global(gindex))
        {
          uint const lindex = dist.get_local(gindex);
          // Check that receive entity global index is shared on current rank
          if (dist.is_shared(lindex))
          {

          }
          else
          {
            error("Global entity %u not shared on rank %u", gindex, pe_rank);
          }
        }
        else
        {
          error("Global entity %u not found on rank %u", gindex, pe_rank);
        }
      }

      recv_buf += recv_count;
      recv_max -= recv_count;
    }
    delete[] sbuf;
  }

  // Check ghost entities adjacency
  {
    Array<uint> * sbuf = new Array<uint> [pe_size];
    uint e_count = 0;
    for (typename E::ghost e(mesh); !e.end(); ++e, ++e_count)
    {
      sbuf[e.owner()].push_back(e.global_index());
    }
    // Check that ghost entities indices were counted correctly
    dolfin_assert(e_count == mesh.topology().num_ghost(edim));

    // Swap local entities to array
    Array<uint> rbuf; rbuf.swap(sbuf[pe_rank]);
    dolfin_assert(sbuf[pe_rank].size() == 0);

    // Exchange entities
    MPI_Status status;
    uint src;
    uint dst;
    uint recv_max;
    for (uint j = 0; j < pe_size; ++j)
    {
      uint s = sbuf[j].size();
      MPI_Reduce(&s, &recv_max, 1, MPI_UNSIGNED, MPI_SUM, j, dist.comm());
      // Check that no duplicate global entity was added
      std::set<uint> global_indices(sbuf[j].begin(), sbuf[j].end());
      if (global_indices.size() != sbuf[j].size())
      {
        error("Duplicate global indices for entities of dimension %u", edim);
      }
    }
    uint const rbuf_size = rbuf.size();
    rbuf.resize(rbuf_size + recv_max);
    uint * recv_buf = rbuf.ptr() + rbuf_size;
    int recv_count;
    _set<uint> recv_idx;
    for (uint j = 1; j < pe_size; ++j)
    {
      src = (pe_rank - j + pe_size) % pe_size;
      dst = (pe_rank + j) % pe_size;

      MPI_Sendrecv(&sbuf[dst][0], sbuf[dst].size(), MPI_UNSIGNED, dst, 0,
                   &recv_buf[0] , recv_max        , MPI_UNSIGNED, src, 0,
                   dist.comm(), &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

      for (int k = 0; k < recv_count; ++k)
      {
        uint const gindex = recv_buf[k];
        // Check that receive entity global index exists on current rank
        if (dist.has_global(gindex))
        {
          uint const lindex = dist.get_local(gindex);
          // Check that receive entity global index is shared on owner rank
          if (!dist.is_shared(lindex))
          {
            error("Global entity %u not shared on rank %u", gindex, pe_rank);
          }
          // Check that sender is listed as adjacent rank (redundant but safer)
          if (dist.get_shared_adj(lindex).count(src) == 0)
          {
            error("Global entity %u not shared with rank %u", gindex, pe_rank);
          }
          // Check that receive entity global index is not ghost on owner rank
          if (dist.is_ghost(lindex))
          {
            error("Global entity %u not ghost on rank %u", gindex, pe_rank);
          }
          else
          {
            recv_idx.insert(gindex);
            // Check count of ghost entity global index (redundant but safer)
            if (recv_idx.count(gindex) > dist.get_shared_adj(lindex).size())
            {
              error("Global ghost entity %u received too many times", gindex);
            }
          }
        }
        else
        {
          error("Global entity %u not found on rank %u", gindex, pe_rank);
        }
      }

      recv_buf += recv_count;
      recv_max -= recv_count;
    }

    delete[] sbuf;
  }
#endif /* HAVE_MPI */
}
//--- TEMPLATE INSTANTIATIONS -------------------------------------------------
template void MPIMeshCommunicator::check<Vertex>(Mesh& Mesh);
template void MPIMeshCommunicator::check<Edge>  (Mesh& Mesh);
template void MPIMeshCommunicator::check<Face>  (Mesh& Mesh);
template void MPIMeshCommunicator::check<Facet> (Mesh& Mesh);
template void MPIMeshCommunicator::check<Cell>  (Mesh& Mesh);
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::check(Mesh& mesh)
{
  uint const tdim = mesh.topology().dim();
  check<Vertex>(mesh);
  if (tdim > 1) check<Edge>(mesh);
  if (tdim > 2) check<Face>(mesh);
  if (tdim > 0) check<Cell>(mesh);
}
//-----------------------------------------------------------------------------
} /* namespace dolfin */
