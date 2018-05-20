// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2016.
// Rewritten to fix distribution bugs and use distributed data class.
//

#include <dolfin/mesh/MeshRenumber.h>

#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/mesh/EntityKey.h>
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
  uint const rank = MPI::rank();
  uint const pe_size = MPI::size();

  /*
   * Renumber vertices: compacting global numbering per rank.
   *
   */

  if (topology.connectivity(0) && !distdata[0].valid_numbering)
  {
    message(1, "MeshRenumber : renumber vertices");
    renumbered = true;
    DistributedData& vdata = distdata[0];
    dolfin_assert(vdata.local_size() == topology.size(0));
    vdata.renumber_global();
    dolfin_assert(vdata.local_size() == topology.size(0));
    vdata.valid_numbering = true;
  }

  /*
   * Renumber edges/faces: determine shared/ghost entities and number.
   *
   */

  for (uint d = 1; d < tdim; ++d)
  {
    if (!(topology.connectivity(d) && !distdata[d].valid_numbering))
    {
      continue;
    }

    message(1, "MeshRenumber : renumber entities of dimension %u", d);
    renumbered = true;
    DistributedData& vdata = distdata[0];
    DistributedData& edata = distdata[d];
    Connectivity const& cve = topology(0, d);
    Connectivity const& cev = topology(d, 0);

    // Set the random seed
    std::srand(MPI::seed());

    // Distributed data size is known: cache arrays are used.
    dolfin_assert(edata.empty());
    dolfin_assert(!edata.is_finalized());
    edata.set_size(topology.size(d));

    //
    uint const num_entity_vertices = cev.max_degree();
    EntityKey key(num_entity_vertices);
    _map<EntityKey, uint> entity_map;
    Array<uint> * sendbuf = new Array<uint> [pe_size];

    // Collect entities with a common adjacent to shared vertices
    _set<uint> adjs;
    bool * used_entities = new bool[topology.size(d)];
    std::fill_n(used_entities, topology.size(d), false);
    for (SharedIterator it(vdata); it.valid(); ++it)
    {
      dolfin_assert(it.index() < cve.order());
      uint const * v_entities = cve(it.index());
      for (uint e = 0; e < cve.degree(it.index()); ++e)
      {
        uint const entity_index = v_entities[e];
        if (used_entities[entity_index] == true)
        {
          continue;
        }
        used_entities[entity_index] = true;

        // Skip entities with a non-shared vertex
        uint const * vertices = cev(entity_index);
        dolfin_assert(vertices[0] != vertices[1]);
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
            vdata.get_global(num_entity_vertices, vertices, key.indices);
            // NOTE: it is important to use set to copy global indices as sort
            //       is called to store indices in increasing order.
            key.set(key.indices, entity_index);
            dolfin_assert(key.indices[0] < key.indices[1]);
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
    delete [] used_entities;

    // Exchange data to mark which entities are shared
    _map<uint,uint> recvmap;

    MPI_Status status;
    uint src;
    uint dst;
    uint sendmax = sendbuf[0].size();
    for (uint j = 1; j < pe_size; ++j)
    {
      sendmax = std::max(sendmax, (uint) sendbuf[j].size());
    }
    uint recvmax = 0;
    MPI::all_reduce<MPI::max>(sendmax, recvmax);
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
        // Data is received in the following order:
        // [ vote, local index, sorted global vertex indices ]
        uint const vote1 = recvbuf[k];
        key.idx = recvbuf[k + 1];
        // Entities are already sorted, just copy them to save sort computation
        std::copy(&recvbuf[k + 2], &recvbuf[k + 2] + num_entity_vertices,
                  key.indices);
        //dolfin_assert(key.indices[0] < key.indices[1]);

        // Beware camembert !
        // Even if the rank is adjacent for all the vertices, the entity may
        // still not be shared: the main bug of the original implementation was
        // that entities with all vertices shared were marked as shared although
        // the condition is only necessary, not sufficient.
        _map<EntityKey, uint>::iterator it = entity_map.find(key);
        if (it != entity_map.end())
        {
          //FIXME: Hash collision possible ?
          uint const local_index = it->first.idx;
          uint const vote0 = it->second;

          // Give ownership to the minimum vote amongst adjacent ranks
          if ((vote1 < vote0) ||
              (vote1 == vote0 && src < edata.get_owner(local_index)))
          {
            // Update vote, map local index and set owner
            it->second = vote1;
            dolfin_assert(it->second == vote1);
            recvmap[local_index] = key.idx;
            dolfin_assert(key.idx == recvbuf[k + 1]);
            edata.set_ghost(local_index, src);
            dolfin_assert(edata.is_ghost(local_index));
            dolfin_assert(edata.get_owner(local_index) == src);
            dolfin_assert(edata.is_shared(local_index));
          }
          edata.set_shared_adj(local_index, src);
          dolfin_assert(edata.is_shared(local_index));
        }
      }
    }

    // Cleanup
    delete[] recvbuf;
    delete[] sendbuf;
    entity_map.clear();

    // Exchange ghost entities
    sendbuf = new Array<uint>[pe_size];
    Array<uint> * ghostid = new Array<uint>[pe_size];
    edata.set_range(cev.order() - edata.num_ghost());
    uint current_index = edata.offset();
    for(uint i = 0; i < cev.order(); ++i)
    {
      if(edata.is_owned(i))
      {
        edata.set_map(i, current_index);
        ++current_index;
      }
      else
      {
        // Enqueue to query global index
        uint const owner = edata.get_owner(i);
        dolfin_assert(owner != rank);
        dolfin_assert(owner < pe_size);
        dolfin_assert(recvmap.count(i) > 0);
        sendbuf[owner].push_back(recvmap.find(i)->second);
        ghostid[owner].push_back(i);
      }
    }

    // At this point the mapping is set for owned entities but not for ghosts
    if (edata.num_shared() < edata.num_ghost())
    {
      error("MeshRenumber : invalid number of entities, shared %u < %u ghost",
            edata.num_shared(), edata.num_ghost());
    }
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
        ///Throw hard error to make sure consistency of data is ensured before
        // renumbering.
        if(!edata.is_shared(recvbuf[k]))
        {
          error("MeshRenumber : received entity %u is not marked as shared",
                recvbuf[k]);
        }
        if(!edata.is_owned(recvbuf[k]))
        {
          error("MeshRenumber : received entity %u is not marked as owned",
                recvbuf[k]);
        }
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
    dolfin_assert(edata.local_size() == topology.size(d));
    edata.valid_numbering = true;
  }

  /*
   * Renumber cells: compacting global numbering per rank.
   *
   */

  if (topology.connectivity(tdim) && !distdata[tdim].valid_numbering)
  {
    message(1, "MeshRenumber : renumber cells");
    renumbered = true;
    DistributedData& cdata = distdata[tdim];
    // Cell numbering is applied automatically at finalized call
    cdata.renumber_global();
    dolfin_assert(cdata.local_size() == topology.size(tdim));
    cdata.valid_numbering = true;
  }

#endif /* HAVE_MPI */

  return renumbered;
}

} /* namespace dolfin */
