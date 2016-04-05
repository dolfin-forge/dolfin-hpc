// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/config/dolfin_config.h>

#include <dolfin/mesh/MeshRenumber.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Face.h>
#include <dolfin/mesh/Vertex.h>

#include <dolfin/common/Array.h>
#include <dolfin/main/MPI.h>

#include <cstdlib>
#include <ctime>
#include <map>
#include <set>
#include <utility>
#include <fstream>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
bool MeshRenumber::renumber(Mesh& mesh)
{
  bool ret = renumber_vertices(mesh);
#ifndef ENABLE_P1_OPTIMIZATIONS
  ret &= renumber_edges(mesh);
  ret &= renumber_faces(mesh);
  ret &= remap_facets(mesh);
#endif
  ret &= renumber_cells(mesh);
  return ret;
}
//-----------------------------------------------------------------------------
#ifdef HAVE_MPI
//-----------------------------------------------------------------------------
bool MeshRenumber::renumber_vertices(Mesh& mesh)
{
  if (mesh.topology().dim() < 1 || !mesh.is_distributed()
      || mesh.distdata()[0].valid_numbering)
  {
    return false;
  }

  DistributedData& distdata0 = mesh.distdata()[0];
//  dolfin_assert(distdata0.is_finalized());

  int const rank = MPI::processNumber();
  int const pe_size = MPI::numProcesses();

  // Update global number of vertices and get offset of local numbering
  uint offset = distdata0.offset();

  DistributedData distdata1;
  for (uint i = 0; i < mesh.numVertices(); ++i)
  {
    if (distdata0.is_owned(i))
    {
      distdata1.set_map(i, offset++);
    }
  }

  // Collect ghosted entities per owner
  Array<uint> *ghost_buff = new Array<uint> [pe_size];
  for (GhostIterator iter(distdata0); !iter.end(); ++iter)
  {
    ghost_buff[iter.owner()].push_back(iter.global_index());
  }

  // Exchange data and set numbering
  MPI_Status status;
  Array<uint> send_buff;
  uint src, dest;
  uint recv_size = distdata0.num_ghost();
  int recv_count, recv_size_gh, send_size;

  for (int i = 0; i < pe_size; ++i)
  {
    send_size = ghost_buff[i].size();
    MPI_Reduce(&send_size, &recv_size_gh, 1, MPI_INT, MPI_SUM, i,
               MPI::DOLFIN_COMM);
  }

  uint *recv_ghost = new uint[recv_size_gh];
  uint *recv_buff = new uint[recv_size];

  for (int j = 1; j < pe_size; ++j)
  {
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&ghost_buff[dest][0], ghost_buff[dest].size(), MPI_UNSIGNED,
                 dest, 1, recv_ghost, recv_size_gh, MPI_UNSIGNED, src, 1,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int k = 0; k < recv_count; ++k)
    {
      uint const local_index = distdata0.get_local(recv_ghost[k]);
      uint const new_global = distdata1.get_global(local_index);
      send_buff.push_back(new_global);
    }

    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, src, 2,
                 recv_buff, recv_size, MPI_UNSIGNED, dest, 2, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int j = 0; j < recv_count; ++j)
    {
      distdata1.set_map(distdata0.get_local(ghost_buff[dest][j]), recv_buff[j]);
    }
    send_buff.clear();
  }

  // Use new numbering
//  distdata0.apply_numbering(0, new_local, new_global);
//  mddata.apply_ownership(0);
//  distdata0[0].finalize();

  delete[] recv_buff;
  delete[] recv_ghost;
  for (int i = 0; i < pe_size; ++i)
  {
    ghost_buff[i].clear();
  }
  delete[] ghost_buff;

  return true;
}
//-----------------------------------------------------------------------------
bool MeshRenumber::renumber_edges(Mesh& mesh)
{
  if (mesh.topology().dim() < 2 || !mesh.is_distributed()
      || mesh.distdata()[1].valid_numbering)
  {
    return false;
  }
  MeshDistributedData& mddata = mesh.distdata();
//  mddata.flush_numbering_data(1);
//  mddata.flush_ownership_data(1);

  int const rank = MPI::processNumber();
  int const pe_size = MPI::numProcesses();

  EdgeKey edgekey;
  std::map<EdgeKey, uint> edge_map;
  std::map<EdgeKey, uint> edge_id;
  std::set<EdgeKey> ghosted_edges;

  Array<uint> send_buff;
  Array<uint> send_buff_id;
  _map<uint,uint> send_mapping;
  _set<uint> used_edge;

  // Collect all potential edges satisfying the necessary condition:
  // all the vertices are shared with a common adjacent.
  // Here we just send edges with all shared vertices.
  srand((uint) time(0) + rank);
  for (SharedIterator sv(mddata[0]); !sv.end(); ++sv)
  {
    Vertex v(mesh, sv.index());
    for (MeshEntityIterator e(v, 1); !e.end(); ++e)
    {
      if (used_edge.count(e->index()) == 0)
      {
        const uint *edge_v = e->entities(0);
        uint w = (sv.index() == edge_v[0] ? edge_v[1] : edge_v[0]);
        // Send only if both vertices are shared
        if (mddata[0].is_shared(w))
        {
          edgekey = edge_key(edge_v[0], edge_v[1]);
          edge_map[edgekey] = e->index();
          edge_id[edgekey] = (uint) rand();
          send_buff.push_back(mddata[0].get_global(edge_v[0]));
          send_buff.push_back(mddata[0].get_global(edge_v[1]));
          send_buff_id.push_back(edge_id[edgekey]);
        }
        used_edge.insert(e->index());
      }
    }
  }

  uint num_ghost = 0;
  // Assign ownership of shared edges
  MPI_Status status;
  uint src, dest;
  int max_un, num_un, max_id, num_id, recv_count;
  num_un = send_buff.size();
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
  num_id = send_buff_id.size();
  MPI_Allreduce(&num_id, &max_id, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
  uint *recv_buff = new uint[max_un];
  uint *recv_buff_id = new uint[max_id];

  for (int j = 1; j < (int) pe_size; ++j)
  {

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff_id[0], num_id, MPI_UNSIGNED, dest, 1, recv_buff_id,
                 max_id, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, recv_buff,
                 max_un, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (uint i = 0; i < (uint) recv_count; i += 2)
    {
      // Check if I have the vertices
      if (mddata[0].has_global(recv_buff[i])
          && mddata[0].has_global(recv_buff[i + 1]))
      {

        // Generate edge key
        edgekey = edge_key(mddata[0].get_local(recv_buff[i]),
                       mddata[0].get_local(recv_buff[i + 1]));

        // Check if I have the corresponding edge
        if (edge_id.count(edgekey))
        {
          // Set entity as ghosted based on voting process and remove from map
          if (recv_buff_id[i >> 1] < edge_id[edgekey]
              || (recv_buff_id[i >> 1] == edge_id[edgekey]
                  && src < rank))
          {
            edge_id.erase(edgekey);
            mddata[1].set_ghost(edge_map[edgekey], src);
            ghosted_edges.insert(edgekey);
            ++num_ghost;
          }
        }
      }
    }
  }
  dolfin_assert(num_ghost == ghosted_edges.size());

  // Update global number of edges and get offset of local numbering
  uint offset = mddata[1].offset();

  send_buff.clear();

  uint num = 0;
  _map<uint,uint> new_local,new_global;
  for (uint i = 0; i < mesh.numEdges(); ++i)
  {
    if (!mddata[1].is_ghost(i))
    {
      new_global[i] = offset++;
      new_local[new_global[i]] = i;
    }
    else
    {
      Edge e(mesh, i);
      const uint *edge_v = e.entities(0);
      send_buff.push_back(mddata[0].get_global(edge_v[0]));
      send_buff.push_back(mddata[0].get_global(edge_v[1]));
      send_mapping[num++] = e.index();
    }
  }
  num_un = send_buff.size();

  //Exchange assigned global numbers
  Array<uint> global_buff;
  for (int j = 1; j < pe_size; ++j)
  {

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, recv_buff,
                 max_un, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int i = 0; i < recv_count; i += 2)
    {
      // Check if I have the vertices
      if (mddata[0].has_global(recv_buff[i])
          && mddata[0].has_global(recv_buff[i + 1]))
      {

        // Generate edge key
        edgekey = edge_key(mddata[0].get_local(recv_buff[i]),
                       mddata[0].get_local(recv_buff[i + 1]));

        // Set entity as shared and add adjacents
        if (edge_id.count(edgekey))
        {
          global_buff.push_back(i >> 1);
          global_buff.push_back(new_global[edge_map[edgekey]]);
          mddata[1].set_shared(edge_map[edgekey]);
          mddata[1].set_shared_adj(edge_map[edgekey], status.MPI_SOURCE);
        }
        else if (ghosted_edges.count(edgekey))
        {
          mddata[1].set_shared_adj(edge_map[edgekey], status.MPI_SOURCE);
        }
      }
    }

    MPI_Sendrecv(&global_buff[0], global_buff.size(), MPI_UNSIGNED, src, 2,
                 recv_buff, max_un, MPI_UNSIGNED, dest, 2, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int i = 0; i < recv_count; i += 2)
    {
      new_global[send_mapping[recv_buff[i]]] = recv_buff[i + 1];
      new_local[recv_buff[i + 1]] = send_mapping[recv_buff[i]];
      mddata[1].set_ghost(send_mapping[recv_buff[i]], status.MPI_SOURCE);
    }
    global_buff.clear();
  }

  // Use new numbering
//  mddata.apply_numbering(1, new_local, new_global);
//  mddata.apply_ownership(1);
  mddata[1].finalize();

  delete[] recv_buff;
  delete[] recv_buff_id;

  return true;
}
//-----------------------------------------------------------------------------
bool MeshRenumber::renumber_faces(Mesh& mesh)
{
  if (mesh.topology().dim() < 3 || !mesh.is_distributed()
      || mesh.distdata()[2].valid_numbering)
  {
    return false;
  }
  MeshDistributedData& mddata = mesh.distdata();
//  mddata.flush_numbering_data(2);
//  mddata.flush_ownership_data(2);

  int const rank = MPI::processNumber();
  int const pe_size = MPI::numProcesses();

  FaceKey facekey;
  std::map<FaceKey, uint> face_map;
  std::map<FaceKey, uint> face_id;
  std::set<FaceKey> ghosted_faces;

  Array<uint> send_buff;
  Array<uint> send_buff_id;
  std::map<uint, uint> send_mapping;

  _set<uint> used_face;

  srand((uint) time(0) + rank);
  BoundaryMesh local_boundary(mesh, BoundaryMesh::interior);
  for (CellIterator bf(local_boundary); !bf.end(); ++bf)
  {
    Face f(mesh, local_boundary.facet_index(*bf));
    send_buffer_face(send_buff, mesh, f);
    facekey = face_key(f);
    face_map[facekey] = f.index();
    face_id[facekey] = (uint) rand();
    send_buff_id.push_back(face_id[facekey]);
    mddata[2].set_shared(f.index());
  }

  uint num_ghost = 0;
  Face f(mesh, 0);
  uint inc = 2 * f.num_entities(0);

  // Assign ownership of shared faces
  MPI_Status status;
  uint src, dest;
  int max_un, num_un, max_id, num_id, recv_count, recv_count_id;
  num_un = send_buff.size();
  MPI_Allreduce(&num_un, &max_un, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
  num_id = send_buff_id.size();
  MPI_Allreduce(&num_id, &max_id, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
  uint *recv_buff = new uint[max_un];
  uint *recv_buff_id = new uint[max_id];
  EdgeKey key;
  for (int j = 1; j < (int) pe_size; ++j)
  {

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff_id[0], num_id, MPI_UNSIGNED, dest, 1, recv_buff_id,
                 max_id, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count_id);

    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, recv_buff,
                 max_un, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    uint ii = 0;
    for (uint i = 0; i < (uint) recv_count; ++ii, i += inc)
    {
      facekey.clear();
      uint num_ok = 0;
      for (uint k = 0; k < inc; k += 2)
      {
        // Check if I have the vertices
        if (mddata[0].has_global(recv_buff[i + k])
            && mddata[0].has_global(recv_buff[i + k + 1]))
        {
          // Generate edge key
          key = edge_key(mddata[0].get_local(recv_buff[i + k]),
                         mddata[0].get_local(recv_buff[i + k + 1]));
          facekey.insert(key);
          ++num_ok;
        }
      }

      if (num_ok < f.num_entities(0))
      {
        continue;
      }

      // Check if I have the corresponding edge
      if (face_id.count(facekey))
      {
        dolfin_assert(face_id.count(facekey));
        if (recv_buff_id[ii] < face_id[facekey]
            || (recv_buff_id[ii] == face_id[facekey] && src < rank))
        {
          face_id.erase(facekey);
          mddata[2].set_ghost(face_map[facekey], src);
          ghosted_faces.insert(facekey);
          ++num_ghost;
        }
      }
    }
  }
  dolfin_assert(num_ghost == ghosted_faces.size());

  // Update global number of faces and get offset of local numbering
  uint offset = mddata[2].offset();

  send_buff.clear();

  uint num = 0;
  _map<uint,uint> new_local,new_global;
  for (uint i = 0; i < mesh.numFaces(); ++i)
  {
    if (!mddata[2].is_ghost(i))
    {
      new_global[i] = offset++;
      new_local[new_global[i]] = i;
    }
    else
    {
      Face f(mesh, i);
      send_mapping[num++] = f.index();
      send_buffer_face(send_buff, mesh, f);
    }
  }

  num_un = send_buff.size();
  //Exchange assigned global numbers
  Array<uint> global_buff;
  for (int j = 1; j < pe_size; ++j)
  {

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff[0], num_un, MPI_UNSIGNED, dest, 1, recv_buff,
                 max_un, MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
    uint ii = 0;
    for (uint i = 0; i < (uint) recv_count; ++ii, i += inc)
    {
      // Check if I have the vertices
      facekey.clear();
      uint num_ok = 0;
      for (uint k = 0; k < inc; k += 2)
      {
        if (mddata[0].has_global(recv_buff[i + k])
            && mddata[0].has_global(recv_buff[i + k + 1]))
        {
          // Generate edge key
          key = edge_key(mddata[0].get_local(recv_buff[i + k]),
                         mddata[0].get_local(recv_buff[i + k + 1]));

          facekey.insert(key);
          ++num_ok;
        }
      }
      if (num_ok < f.num_entities(0))
      {
        continue;
      }

      if (face_id.count(facekey))
      {
        global_buff.push_back(ii);
        global_buff.push_back(new_global[face_map[facekey]]);
        mddata[2].set_shared(face_map[facekey]);
        mddata[2].set_shared_adj(face_map[facekey], src);
      }
      else if (ghosted_faces.count(facekey))
      {
        mddata[2].set_shared_adj(face_map[facekey], status.MPI_SOURCE);
      }
    }

    MPI_Sendrecv(&global_buff[0], global_buff.size(), MPI_UNSIGNED, src, 2,
                 recv_buff, max_un, MPI_UNSIGNED, dest, 2, MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int i = 0; i < recv_count; i += 2)
    {
      new_global[send_mapping[recv_buff[i]]] = recv_buff[i + 1];
      new_local[recv_buff[i + 1]] = send_mapping[recv_buff[i]];
      mddata[2].set_ghost(send_mapping[recv_buff[i]], status.MPI_SOURCE);
    }
    global_buff.clear();
  }

  // Use new numbering
//  mddata.apply_numbering(2, new_local, new_global);
//  mddata.apply_ownership(2);
  mddata[2].finalize();

  delete[] recv_buff;
  delete[] recv_buff_id;

  return true;
}
//-----------------------------------------------------------------------------
bool MeshRenumber::renumber_cells(Mesh& mesh)
{
  uint const tdim = mesh.topology().dim();
  if (!mesh.is_distributed() || mesh.distdata()[tdim].valid_numbering)
  {
    return false;
  }

  MeshDistributedData& mddata = mesh.distdata();
  uint offset = mddata[tdim].offset();

  _map<uint,uint> new_local;
  _map<uint,uint> new_global;
  for (uint i = 0; i < mesh.numCells(); ++i)
  {
    new_global[i] = offset++;
    new_local[new_global[i]] = i;
  }

  // Use new numbering
//  mddata.apply_numbering(tdim, new_local, new_global);
//  mddata.apply_ownership(tdim);
  mddata[tdim].finalize();

  return true;
}
//-----------------------------------------------------------------------------
bool MeshRenumber::remap_facets(Mesh& mesh)
{
  uint const facetdim = mesh.topology().dim() - 1;
  if (!mesh.is_distributed())// || mesh.distdata().has_valid_mapping(facetdim))
  {
    return false;
  }
  remap_shared_entities(mesh, facetdim);
  return true;
}
//-----------------------------------------------------------------------------
std::pair<uint, uint> MeshRenumber::edge_key(uint id1, uint id2)
{
  if (id2 < id1)
  {
    EdgeKey key(id2, id1);
    return key;
  }
  else
  {
    EdgeKey key(id1, id2);
    return key;
  }

}
//-----------------------------------------------------------------------------
std::set<std::pair<uint, uint> > MeshRenumber::face_key(Face& f)
{
  const uint *face_v = f.entities(0);
  FaceKey fk;
  fk.insert(edge_key(face_v[0], face_v[1]));
  fk.insert(edge_key(face_v[1], face_v[2]));

  switch (f.num_entities(0))
    {
    case 3:
      fk.insert(edge_key(face_v[2], face_v[0]));
      break;
    case 4:
      fk.insert(edge_key(face_v[2], face_v[3]));
      fk.insert(edge_key(face_v[3], face_v[0]));
      break;
    default:
      error("Unkown entity");
      break;
    }

  return fk;
}
//-----------------------------------------------------------------------------
void MeshRenumber::send_buffer_face(Array<uint>& send_buff, Mesh& mesh, Face& f)
{
  MeshDistributedData& mddata = mesh.distdata();
  const uint *face_v = f.entities(0);
  send_buff.push_back(mddata[0].get_global(face_v[0]));
  send_buff.push_back(mddata[0].get_global(face_v[1]));

  send_buff.push_back(mddata[0].get_global(face_v[1]));
  send_buff.push_back(mddata[0].get_global(face_v[2]));

  switch (f.num_entities(0))
    {
    case 3:
      send_buff.push_back(mddata[0].get_global(face_v[2]));
      send_buff.push_back(mddata[0].get_global(face_v[0]));
      break;
    case 4:
      send_buff.push_back(mddata[0].get_global(face_v[2]));
      send_buff.push_back(mddata[0].get_global(face_v[3]));

      send_buff.push_back(mddata[0].get_global(face_v[3]));
      send_buff.push_back(mddata[0].get_global(face_v[0]));
      break;
    default:
      error("Unkown entity");
      break;
    }
}
//-----------------------------------------------------------------------------
void MeshRenumber::remap_shared_entities(Mesh& mesh, uint const dim)
{
  /*
  if(!mesh.is_distributed())
  {
    return;
  }
  MeshDistributedData& distdata = mesh.distdata();
//  mesh.distdata().flush_mapping(dim);

  //
  uint rank = dolfin::MPI::processNumber();
  uint pe_size = dolfin::MPI::numProcesses();

  // Array of global indices of shared entities to be sent
  Array<uint> * sendbuf = new Array<uint> [pe_size];

  // Mapping from global entity indices to shared iterator indices
  uint * idx = new uint[pe_size];

  // Group shared entities by adjacent rank
  _map<uint,uint> * shared_idx = new _map<uint,uint> [pe_size];
  std::fill_n(idx, pe_size, 0);
  for (SharedIterator sh(distdata[dim]); !sh.end(); ++sh)
  {
    uint const glb = distdata[dim].get_global(sh.index());
    _set<uint> const& adj = distdata[dim].get_shared_adj(sh.index());
    for (_set<uint>::const_iterator it = adj.begin(); it != adj.end(); ++it)
    {
      sendbuf[*it].push_back(glb);
      shared_idx[*it].insert(std::pair<uint,uint>(glb, idx[*it]++));
    }
  }

  // Group ghost entities by owner
  _map<uint,uint> * ghost_idx = new _map<uint,uint> [pe_size];
  std::fill_n(idx, pe_size, 0);
  for (GhostIterator gh(distdata[dim]); !gh.end(); ++gh)
  {
    uint const glb = distdata[dim].get_global(gh.index());
    uint owner = distdata[dim].get_owner(gh.index());
    ghost_idx[owner].insert(std::pair<uint, uint>(glb, idx[owner]++));
  }

  //
  delete[] idx;

  //
  MPI_Status status;
  int src = 0;
  int dest = 0;
  int num_send_entities = 0;
  int num_recv_entities = 0;
  int recv_count = 0;

  // Exchange global indices to know ordering of adjacent rank
  // Set recv buffer size to number of shared entities to avoid reallocation
  uint * recvbuf = new uint[mesh.distdata()[dim].num_shared()];
  for (int j = 1; j < (int) pe_size; ++j)
  {
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    num_send_entities = sendbuf[dest].size();
    num_recv_entities = sendbuf[src].size();

    MPI_Sendrecv(&sendbuf[dest][0], num_send_entities, MPI_UNSIGNED, dest, 1,
                 &recvbuf[0], num_recv_entities, MPI_UNSIGNED, src, 1,
                 dolfin::MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    if (recv_count != num_recv_entities)
    {
      error("Mismatch between the number of sent and received entities.");
    }

    // Skip non-adjacent ranks
    if (distdata[dim].get_adj_ranks().count(src) == 0) continue;

    // Shared and ghost mappings
    Array<uint>& shared_mapping0 = distdata.get_shared_mapping_to(src, dim);
    shared_mapping0.resize(shared_idx[src].size());
    Array<uint>& shared_mapping1 = distdata.get_shared_mapping_from(src, dim);
    shared_mapping1.resize(shared_idx[src].size());
    Array<uint>& ghost_mapping0 = distdata.get_ghost_mapping_to(src, dim);
    ghost_mapping0.resize(ghost_idx[src].size());
    Array<uint>& ghost_mapping1 = distdata.get_ghost_mapping_from(src, dim);
    ghost_mapping1.resize(ghost_idx[src].size());
    uint ghost_entity = 0;
    for (uint entity = 0; entity < (uint) num_recv_entities; ++entity)
    {
      uint glb = recvbuf[entity];
      //
      dolfin_assert(distdata[dim].has_global(glb));
      _map<uint, uint>::const_iterator sit = shared_idx[src].find(glb);
      dolfin_assert(sit != shared_idx[src].end());

      // Fill shared_mappings between local and adjacent shared ordering
      // Anytime facets are send/received from a given rank, they are ordered
      // according to the given arrays.
      shared_mapping0[sit->second] = entity;  // send (local to adjacent)
      shared_mapping1[entity] = sit->second;  // recv (adjcent to local)

      _map<uint, uint>::const_iterator git = ghost_idx[src].find(glb);
      if (git != ghost_idx[src].end())
      {
        //
        dolfin_assert(
            distdata[dim].is_ghost(distdata[dim].get_local(glb))
                && (distdata[dim].get_owner(distdata[dim].get_local(glb))
                    == (uint ) src));

        // Fill the ghost mapping between local ghost iterator ordering and the
        // ordering of corresponding *shared* entities for the owner rank.
        ghost_mapping0[git->second] = ghost_entity;  // send
        ghost_mapping1[ghost_entity] = git->second;  // recv
        ++ghost_entity;
      }

      //
      dolfin_assert(shared_mapping0.size() == shared_mapping1.size());
      dolfin_assert(ghost_mapping0.size() == ghost_mapping1.size());

    }

  }
  delete[] recvbuf;
  delete[] ghost_idx;
  delete[] shared_idx;
  delete[] sendbuf;
  */
}
//-----------------------------------------------------------------------------
#else
//-----------------------------------------------------------------------------
bool MeshRenumber::renumber_vertices(Mesh& mesh)
{
  return true;
}
//-----------------------------------------------------------------------------
bool MeshRenumber::renumber_edges(Mesh& mesh)
{
  return true;
}
//-----------------------------------------------------------------------------
bool MeshRenumber::renumber_faces(Mesh& mesh)
{
  return true;
}
//-----------------------------------------------------------------------------
bool MeshRenumber::renumber_cells(Mesh& mesh)
{
  return true;
}
//-----------------------------------------------------------------------------
bool MeshRenumber::remap_facets(Mesh& mesh)
{
  return true;
}
//-----------------------------------------------------------------------------
#endif

}

