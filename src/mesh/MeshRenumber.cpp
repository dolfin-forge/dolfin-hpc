// Copyright (C) 2008 Niclas Jansson.
// Copyright (C) 2016 Aurelien Larcher (rewrite).
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/mesh/MeshRenumber.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/mesh/EntityKey.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
bool MeshRenumber::renumber(MeshTopology& topology)
{
  if (!topology.is_distributed())
  {
    return false;
  }

  bool renumbered = false;

#ifdef HAVE_MPI

  MeshDistributedData& distdata = topology.distdata();

  uint const tdim = topology.dim();
  uint const rank = MPI::processNumber();
  uint const pe_size = MPI::numProcesses();

  /*
   * Renumber vertices: compacting global numbering per rank.
   *
   */

  if (topology.entities_exist(0) && !distdata[0].valid_numbering)
  {
    DistributedData& vdata = distdata[0];
    vdata.finalize();
    vdata.renumber_global();
    vdata.valid_numbering = true;
  }

  /*
   * Renumber edges/faces: determine shared/ghost entities and number.
   *
   */

  for (uint d = 0; d < tdim; ++d)
  {
    if (!topology.entities_exist(d) || distdata[d].valid_numbering)
    {
      continue;
    }

    DistributedData& vdata = distdata[0];
    DistributedData& edata = distdata[d];
    edata.clear();
    MeshConnectivity const& cve = topology(0, d);
    MeshConnectivity const& cev = topology(d, 0);

    //
    uint const num_entity_vertices = cev.max_connections();
    EntityKey key(num_entity_vertices);
    _map<EntityKey, uint> entity_map;
    Array<uint> * sendbuf = new Array<uint> [pe_size];

    // Collect entities with a common adjacent to shared vertices
    _set<uint> adjs;
    _set<uint> used_entities;
    for (SharedIterator it(vdata); !it.end(); ++it)
    {
      dolfin_assert(it.index() < cve.num_entities());
      uint const * v_entities = cve(it.index());
      for (uint e = 0; e < cve.size(it.index()); ++e)
      {
        uint const entity_index = v_entities[e];

        // Skip entities with a non-shared vertex
        uint const * vertices = cev(entity_index);
        bool all_shared = true;
        for (uint v = 0; v < num_entity_vertices; ++v)
        {
          if(it.index() != vertices[v] && !vdata.is_shared(vertices[v]))
          {
            all_shared = false;
            break;
          }
        }

        // Append to send buffer for each adjacent
        if(all_shared)
        {
          //FIXME: randomness may be harmful
          uint const vote = std::rand();
          vdata.get_common_adj(num_entity_vertices, vertices, adjs);
          if (adjs.size() > 0)
          {
            key.idx = entity_index;
            for (uint v = 0; v < num_entity_vertices; ++v)
            {
              key.indices[v] = vdata.get_global(vertices[v]);
            }
            entity_map[key] = vote;
            for (_set<uint>::const_iterator a = adjs.begin(); a != adjs.end();
                 ++a)
            {
              sendbuf[*a].push_back(vote);
              sendbuf[*a].push_back(entity_index);
              for (uint v = 0; v < num_entity_vertices; ++v)
              {
                sendbuf[*a].push_back(key.indices[v]);
              }
            }
          }
        }
      }
    }
    used_entities.clear();

    // Exchange data to mark which entities are shared
    _set<uint> shared;
    _map<uint,uint> recvmap;

    MPI_Status status;
    uint src;
    uint dst;
    uint sendmax = 0;
    for (uint j = 0; j < pe_size; ++j)
    {
      sendmax = std::max(sendmax, (uint) sendbuf[j].size());
    }
    uint recvmax = 0;
    MPI::numGlobalMax(sendmax, recvmax);
    uint * recvbuf = new uint[recvmax];
    int recvcount;
    for (uint j = 1; j < pe_size; ++j)
    {
      src = (rank - j + pe_size) % pe_size;
      dst = (rank + j) % pe_size;

      MPI_Sendrecv(&sendbuf[dst][0], sendbuf[dst].size(), MPI_UNSIGNED, dst, 1,
                   &recvbuf[0], recvmax, MPI_UNSIGNED, src, 1,
                   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recvcount);

      for (int k = 0; k < recvcount; k+=(2 + num_entity_vertices))
      {
        uint const vote1 = recvbuf[k];
        key.set(&recvbuf[k + 2], recvbuf[k + 1]);
        _map<EntityKey, uint>::iterator it = entity_map.find(key);
        // Even if the rank is adjacent for all the vertices, the entity may
        // still not be shared, beware camembert !
        if (it != entity_map.end())
        {
          uint const local_index = it->first.idx;
          uint const vote0 = it->second;
          // Give ownership to the minimum vote amongst adjacent ranks
          if ((vote1 < vote0) || (vote1 == vote0 && src < rank))
          {
            // Update vote, map local index and set owner
            it->second = vote1;
            recvmap[local_index] = key.idx;
            dolfin_assert(key.idx == recvbuf[k + 1]);
            edata.set_ghost(local_index, src);
          }
          else
          {
            edata.set_shared_adj(local_index, src);
          }
          shared.insert(local_index);
        }
      }
    }
    dolfin_assert(edata.capacity() == cev.num_entities());

    // Cleanup
    delete[] recvbuf;
    delete[] sendbuf;
    entity_map.clear();

    // Exchange ghost entities
    sendbuf = new Array<uint>[pe_size];
    Array<uint> * ghostid = new Array<uint>[pe_size];
    uint const num_owned = cev.num_entities() - edata.num_ghost();
    uint offset = 0;
    MPI::processOffset(num_owned, offset);
    uint current_index = 0;
    for(uint i = 0; i < cev.num_entities(); ++i)
    {
      if(edata.is_owned(i))
      {
        edata.set_map(i, offset + current_index);
        ++current_index;
      }
      else
      {
        // Enqueue to query global index
        uint const owner = edata.get_owner(i);
        dolfin_assert(owner != rank);
        dolfin_assert(recvmap.count(i) > 0);
        sendbuf[owner].push_back(recvmap[i]);
        ghostid[owner].push_back(i);
      }
    }

    // At this point the mapping is set for owned entities but not for ghosts
    recvmax = edata.num_shared() - edata.num_ghost();
    recvbuf = (recvmax == 0 ? NULL : new uint[recvmax]);
    uint * sendbuf_back = (recvmax == 0 ? NULL : new uint[recvmax]);
    uint const num_ghosts = edata.num_ghost();
    uint * recvbuf_back = (num_ghosts == 0 ? NULL : new uint[num_ghosts]);
    for (uint j = 1; j < pe_size; ++j)
    {
      src = (rank - j + pe_size) % pe_size;
      dst = (rank + j) % pe_size;

      MPI_Sendrecv(&sendbuf[dst][0], sendbuf[dst].size(), MPI_UNSIGNED, dst, 1,
                   &recvbuf[0], recvmax, MPI_UNSIGNED, src, 1,
                   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recvcount);

      for (int k = 0; k < recvcount; ++k)
      {
        dolfin_assert(shared.count(recvbuf[k]));
        dolfin_assert(edata.is_owned(recvbuf[k]));
        dolfin_assert(edata.is_shared(recvbuf[k]));
        dolfin_assert(edata.get_shared_adj(recvbuf[k]).count(src) > 0);
        sendbuf_back[k] = edata.get_global(recvbuf[k]);
      }

      MPI_Sendrecv(&sendbuf_back[0], recvcount, MPI_UNSIGNED, src, 2,
                   &recvbuf_back[0], sendbuf[dst].size(), MPI_UNSIGNED, dst, 2,
                   MPI::DOLFIN_COMM, &status);

      for (int k = 0; k < (int) sendbuf[dst].size(); ++k)
      {
        edata.set_map(ghostid[dst][k], recvbuf_back[k]);
      }
    }

    delete[] recvbuf_back;
    delete[] sendbuf_back;
    delete[] recvbuf;
    delete[] ghostid;
    delete[] sendbuf;

    //
    edata.finalize();
    edata.valid_numbering = true;
  }

  /*
   * Renumber cells: compacting global numbering per rank.
   *
   */

  if (topology.entities_exist(tdim) && !distdata[tdim].valid_numbering)
  {
    DistributedData& cdata = distdata[tdim];
    cdata.finalize();
    cdata.renumber_global();
    cdata.valid_numbering = true;
  }

#endif /* HAVE_MPI */

  return renumbered;
}

} /* namespace dolfin */
