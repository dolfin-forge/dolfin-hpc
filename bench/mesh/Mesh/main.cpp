#include <dolfin.h>

#include "MeshChecks.h"

using namespace dolfin;

int main(int argc, char** argv)
{
  int ret = 0;

  //---------------------------------------------------------------------------
  Test t(argc, argv);
  Mesh mesh(t.args.mesh_file);
  bool throw_error = true;

  t.begin("Check distributed data");
  {
    if (mesh.distdata().check(throw_error))
    {
      error("Invalid distributed data");
    }
  }
  t.end();

  t.begin("Check common adjacents");
  {
    // Every shared entity is defined such that lower dimensional entities
    // have a common adjacent
    uint const tdim = mesh.topology().dim();
    for (uint e0dim = 1; e0dim <= tdim; ++e0dim)
    {
      message("Check shared entities adjacents for dim %d", e0dim);
      MeshDistributedData& distdata = mesh.distdata();
      uint const e1dim = e0dim - 1;
      for (SharedIterator se(distdata[e0dim]); !se.end(); ++se)
      {
        MeshEntity e0(mesh, e0dim, se.index());
        MeshEntityIterator e1(e0, e1dim);
        _set<uint> cmn_adjs = distdata[e1dim].get_shared_adj(e1->index());
        if (cmn_adjs.empty())
        {
          error("First sub-entity of shared entity has zero adjacents");
        }
        for (++e1; !e1.end(); ++e1)
        {
          _set<uint> const& ae1 = distdata[e1dim].get_shared_adj(e1->index());
          for(_set<uint>::const_iterator ai = cmn_adjs.begin();
              ai != cmn_adjs.end(); ++ai )
          {
            if(ae1.count(*ai) == 0)
            {
              cmn_adjs.erase(ai);
            }
          }
        }
        if (cmn_adjs.empty())
        {
          error("Shared entity does not have common adjacents");
        }
        else
        {
          _set<uint> const& ae0 = distdata.get_shared_adj(e0);
          if(ae0.size() > cmn_adjs.size())
          {
            error("Number of adjacents of shared entity (%u, %u) greater than "
                  "the number of common adjacents of the lower dimensional "
                  "entities.", e0.index(), e0dim);
          }
          else
          {
            for(_set<uint>::const_iterator ai = ae0.begin(); ai != ae0.end();
                ++ai )
            {
              if(cmn_adjs.count(*ai) == 0)
              {
                error("Adjacent ranks set of shared entity (%u, %u) is not a "
                      "subset of common adjacents of the lower dimensional "
                      "entities.", e0.index(), e0dim);
              }
            }
          }
        }
      }
    }
  }
  t.end();

  t.begin("Check entities distribution");
  {
    uint const tdim = mesh.topology().dim();
    for (uint edim = 0; edim < tdim; ++edim)
    {
      message("Check shared and ghost entities distribution for dim %d", edim);
      uint rank = dolfin::MPI::rank();
      uint pe_size = dolfin::MPI::size();
      MeshDistributedData& distdata = mesh.distdata();

      // Shared
      uint const num_shared = distdata[edim].num_shared();
      Array<uint> * sharedbuf = new Array<uint> [pe_size];
      for (SharedIterator se(distdata[edim]); !se.end(); ++se)
      {
        _set<uint> const& adjs = se.adj();
        uint glb_id = distdata[edim].get_global(se.index());
        for(_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it)
        {
          sharedbuf[*it].push_back(glb_id);
          dolfin_assert(sharedbuf[*it].size() <= num_shared);
        }
      }

      // Ghosts
      uint const num_ghost = distdata[edim].num_ghost();
      uint const num_shared_owned = distdata[edim].num_shared() - num_ghost;
      Array<uint> * ghostdbuf = new Array<uint> [pe_size];
      for (GhostIterator ge(distdata[edim]); !ge.end(); ++ge)
      {
        uint const owner = ge.owner();
        uint glb_id = distdata[edim].get_global(ge.index());
        ghostdbuf[owner].push_back(glb_id);
        dolfin_assert(ghostdbuf[owner].size() <= num_ghost);
      }

      //
      MPI_Status status;
      int src = 0;
      int dest = 0;

      int srecv_count = num_shared;
      int smaxrev_count = 0;
      MPI_Allreduce(&srecv_count, &smaxrev_count, 1, MPI_INT, MPI_MAX,
                    dolfin::MPI::DOLFIN_COMM);
      dolfin_assert(smaxrev_count > 0);
      uint * srecv_buff = new uint[smaxrev_count];
      _set<uint> shared_recv;
      int grecv_count = num_ghost;
      int gmaxrev_count = 0;
      _set<uint> owned_recv;
      MPI_Allreduce(&grecv_count, &gmaxrev_count, 1, MPI_INT, MPI_MAX,
                    dolfin::MPI::DOLFIN_COMM);
      dolfin_assert(gmaxrev_count > 0);
      uint * grecv_buff = new uint[gmaxrev_count];
      for (uint p = 1; p < pe_size; ++p)
      {
        src = (rank - p + pe_size) % pe_size;
        dest = (rank + p) % pe_size;

        //
        MPI_Sendrecv(&sharedbuf[dest][0], sharedbuf[dest].size(), MPI_UNSIGNED,
                     dest, 0, srecv_buff, smaxrev_count, MPI_UNSIGNED, src, 0,
                     dolfin::MPI::DOLFIN_COMM, &status);
        MPI_Get_count(&status, MPI_UNSIGNED, &srecv_count);

        uint num_shared_with = 0;
        message("Check shared entities from %d", src);
        for (uint ii = 0; ii < srecv_count; ++ii)
        {
          uint glb_id = srecv_buff[ii];
          // Check if rank has the entity
          if (!distdata[edim].has_global(glb_id))
          {
            error("Unknown shared entity %u received from %d", glb_id, src);
          }
          shared_recv.insert(glb_id);
          // Check if the entity is shared
          uint loc_id = distdata[edim].get_local(glb_id);
          if (!distdata[edim].is_shared(loc_id))
          {
            error("Entity %u is not marked as shared", glb_id);
          }
          // Check if the rank is listed as adjacent
          if (distdata[edim].get_shared_adj(loc_id).count(src) == 0)
          {
            if (distdata[edim].is_ghost(loc_id))
            {
              error("Ghosted entity %u is not shared with %u", glb_id, src);
            }
            else
            {
              error("Owned entity %u is not shared with %u", glb_id, src);
            }
          }
          else
          {
            ++num_shared_with;
          }
        }
        if (distdata.num_shared_with(src, edim) > 0)
        {
          if (num_shared_with != distdata.num_shared_with(src, edim))
          {
            error("Inconsistent number of entities shared with %d: "
                  "(count) %u != %u (dist)",
                  src, num_shared_with, distdata.num_shared_with(src, edim));
          }
        }
        else
        {
          warning("Number of entity shared with %d is not computed", src);
        }

        //
        MPI_Sendrecv(&ghostdbuf[dest][0], ghostdbuf[dest].size(), MPI_UNSIGNED,
                     dest, 1, grecv_buff, gmaxrev_count, MPI_UNSIGNED, src, 1,
                     dolfin::MPI::DOLFIN_COMM, &status);
        MPI_Get_count(&status, MPI_UNSIGNED, &grecv_count);

        message("Check ghost entities from %d", src);
        for (uint ii = 0; ii < grecv_count; ++ii)
        {
          uint glb_id = grecv_buff[ii];
          // Check if it has the vertex
          if (!distdata[edim].has_global(glb_id))
          {
            error("Unknown owned entity %u", glb_id);
          }
          owned_recv.insert(glb_id);
          uint loc_id = distdata[edim].get_local(glb_id);
          if (distdata[edim].is_ghost(loc_id))
          {
            error("Entity %u is not marked as owned", glb_id);
          }
          if (distdata[edim].get_shared_adj(loc_id).count(src) == 0)
          {
            error("Owned entity %u is not shared with %u", glb_id, src);
          }
        }
      }
      if (owned_recv.size() != num_shared_owned)
      {
        error("Mismatch between the number of entities owned and received: "
              "%u != %u", num_shared_owned, owned_recv.size());
      }
      delete[] grecv_buff;
      delete[] srecv_buff;
      delete[] ghostdbuf;
      delete[] sharedbuf;
    }
  }

  t.begin("Check mesh entities sharedness");
  {
    //
    uint const tdim = mesh.topology().dim();
    MeshDistributedData& distdata = mesh.distdata();
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      bool csh = c->is_shared();
      for (uint e1dim = 1; e1dim <= tdim; ++e1dim)
      {
        mesh.init(c->dim(), e1dim);
        for (MeshEntityIterator e1(*c, e1dim); !e1.end(); ++e1)
        {
          if (e1->is_shared())
          {
            _set<uint> e1a = distdata.get_shared_adj(*e1);
            if(e1a.empty())
            {
              error("Adjacents are empty for entity (%u, %u)", e1->index(),
                  e1->dim());
            }

            uint num_shared = 0;
            mesh.init(e1dim, e1dim - 1);
            for (MeshEntityIterator e0(*e1, e1->dim() - 1); !e0.end(); ++e0)
            {
              if (e0->is_shared())
              {
                ++num_shared;
                _set<uint> const& e0a = distdata.get_shared_adj(*e0);

                // All e0 entities have a superset of the e1 adjacent set
                if (e0a.empty())
                {
                  error("Entity (%u, %u) adjacent set is empty.", e0->index(),
                      e0->dim());
                }
                if (e0a.size() < e1a.size())
                {
                  error("Entity (%u, %u) adjacent set smaller than entity "
                      "(%u, %u) set: %u < %u", e0->index(), e0->dim(),
                      e1->index(), e1->dim(), e0a.size(), e1a.size());
                }
                for (_set<uint>::iterator it = e1a.begin(); it != e1a.end();
                    ++it)
                {
                  if(e0a.count(*it) == 0)
                  {
                    dolfin::LogManager::logger().verbose(0);
                    message("Entity (%u, %u) does not have adjacent %u from entity "
                        "(%u, %u)", e0->index(), e0->dim(), *it, e1->index(),
                        e1->dim());
                    std::stringstream ss1;
                    for (_set<uint>::const_iterator it1 = e1a.begin();
                        it1 != e1a.end(); ++it1)
                    {
                      ss1 << *it1 << " ";
                    }
                    message("e1 : is_ghost = %u", e1->is_ghost());
                    message("e1 : %s", ss1.str().c_str());
                    std::stringstream ss0;
                    for (_set<uint>::const_iterator it0 = e0a.begin();
                        it0 != e0a.end(); ++it0)
                    {
                      ss0 << *it0 << " ";
                    }
                    message("e0 : is_ghost = %u", e0->is_ghost());
                    message("e0 : %s", ss0.str().c_str());
                    dolfin::LogManager::logger().silence();
                  }
                }
              }
            }
            // If an entity is shared then all the lower dimensional entities are
            ret &= (num_shared == e1->num_entities(e1->dim() - 1));
            // ! An entity if shared only if they share a common adjacent !
          }
        }
      }
    }
  }
  t.end();

  t.begin("Check interior boundary entities");
  {
    uint const tdim = mesh.topology().dim();
    uint failed = 0;
    for (uint i = 0; i < tdim; ++i)
    {
      bool ok = true;
      ok &= interior_boundary_entities_check(mesh, i, throw_error);
      if (!ok)
      {
        ++failed;
      }
    }
    message("Failed : %d", failed);
    ret += failed;
  }
  t.end();

  t.begin("Check exterior boundary entities");
  {
    uint const tdim = mesh.topology().dim();
    uint failed = 0;
    for (uint i = 0; i < tdim; ++i)
    {
      bool ok = true;
      ok &= exterior_boundary_entities_check(mesh, i, throw_error);
      if (!ok)
      {
        ++failed;
      }
    }
    message("Failed : %d", failed);
    ret += failed;
  }
  t.end();

  return ret;
}

