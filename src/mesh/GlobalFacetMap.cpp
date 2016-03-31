// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-02-25
// Last changed: 2009-03-03

#include <dolfin/mesh/GlobalFacetMap.h>

#include <dolfin/main/MPI.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/EntityKey.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
GlobalFacetMap::GlobalFacetMap(Mesh& mesh) :
    mesh_(mesh),
    tdim_(mesh_.topology().dim())
{
  init();
}
//-----------------------------------------------------------------------------
GlobalFacetMap::~GlobalFacetMap()
{
}
//-----------------------------------------------------------------------------
void GlobalFacetMap::init()
{
  if (!mesh_.is_distributed())
  {
    return;
  }

  message(1, "GlobalFacetMap: compute facet map");
  MeshDistributedData const& distdata = mesh_.distdata();
  shared_facets_.clear();

  if (tdim_ == 1)
  {
    for (MeshSharedIterator it(distdata, 0); !it.end(); ++it)
    {
      shared_facets_.insert(it.index());
    }
  }
  else
  {
    // Iterate over all facets connected to the shared vertices to collect
    // facets on the interprocess boundary and their neighbours
    for (MeshSharedIterator it(distdata, 0); !it.end(); ++it)
    {
      Vertex v(mesh_, it.index());
      for (FacetIterator f(v); !f.end(); ++f)
      {
        if (f->num_entities(tdim_) == 1)
        {
          shared_facets_.insert(f->index());
        }
      }
    }

#ifdef HAVE_MPI
    // Collect facets for which all vertices are shared
    uint const pe_size = MPI::numProcesses();
    uint const rank = MPI::processNumber();
    Array<uint> * sendbuf_facets = new Array<uint> [pe_size];
    Array<uint> * sendbuf_vertices = new Array<uint> [pe_size];
    //FIXME: Cannot work with heterogeneous mesh since the data packing is not
    //       constant, maybe use the maximum
    uint const num_facet_vertices = mesh_.type().num_vertices(tdim_ - 1);
    uint * facet_vertices = new uint[num_facet_vertices];
    _set<uint> adj;
    for (_set<uint>::iterator it = shared_facets_.begin();
         it != shared_facets_.end(); ++it)
    {
      Facet f(mesh_, *it);

      if (f.has_all_vertices_shared())
      {
        adj = distdata.get_shared_adj(f.entities(0)[0], 0);
        for (uint v = 1; v < f.num_entities(0); ++v)
        {
          _set<uint> const& adjx = distdata.get_shared_adj(f.entities(0)[v], 0);
          for(_set<uint>::iterator it = adj.begin(); it != adj.end();)
          {
            if(adjx.count(*it) == 0)
            {
              adj.erase(it++);
            }
            else
            {
              ++it;
            }
          }
        }
        for(_set<uint>::const_iterator it = adj.begin(); it != adj.end();
            ++it)
        {
          for (uint v = 0; v < f.num_entities(0); ++v)
          {
            sendbuf_vertices[*it].push_back(distdata.get_global(f.entities(0)[v], 0));
          }
          sendbuf_facets[*it].push_back(f.index());
        }
      }
    }
    // Clear facets and re-add only if a corresponding facet has been found
    shared_facets_.clear();

    //
    MPI_Status status;
    uint src;
    uint dest;
    uint max_sendsize = 0;
    uint max_rvalsize = 0;
    for (uint j = 0; j < pe_size; ++j)
    {
      max_sendsize = std::max(max_sendsize, (uint) sendbuf_vertices[j].size());
      max_rvalsize = std::max(max_rvalsize, (uint) sendbuf_facets[j].size());
    }
    uint max_recvcount;
    MPI::numGlobalMax(max_sendsize, max_recvcount);
    uint *recvbuf_facets = new uint[max_recvcount];
    uint max_rvalcount;
    MPI::numGlobalMax(max_rvalsize, max_rvalcount);
    uint *sendbuf_rval = new uint[max_rvalcount];
    uint *recvbuf_rval = new uint[max_rvalcount];
    int recv_count = 0;
    for (uint j = 1; j < pe_size; ++j)
    {
      src = (rank - j + pe_size) % pe_size;
      dest = (rank + j) % pe_size;

      MPI_Sendrecv(&sendbuf_vertices[dest][0], sendbuf_vertices[dest].size(),
                   MPI_UNSIGNED, dest, 1, &recvbuf_facets[0], max_recvcount,
                   MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

      uint kfacet = 0;
      uint num_shared_facets = 0;
      for (int k = 0; k < recv_count; k += num_facet_vertices, ++kfacet)
      {
        // If the target rank owns all the vertices then
        // check if there exists a facet composed of these vertices.
        for (uint v = 0; v < num_facet_vertices; ++v)
        {
          dolfin_assert(distdata.has_global(recvbuf_facets[k + v], 0));
          facet_vertices[v] = distdata.get_local(recvbuf_facets[k + v], 0);
        }
        // Pick first vertex and check if one facet matches
        Vertex v0(mesh_, facet_vertices[0]);
        for (FacetIterator f(v0); !f.end(); ++f)
        {
          if (f->num_entities(tdim_) != 1)
          {
            continue;
          }
          uint matching_vertices = 0;
          for (VertexIterator fv(*f); !fv.end(); ++fv)
          {
            for (uint v = 0; v < num_facet_vertices; ++v)
            {
              if (fv->index() == facet_vertices[v])
              {
                ++matching_vertices;
              }
            }
          }
          // Check if a matching facet is found
          if (matching_vertices == num_facet_vertices)
          {
            sendbuf_rval[num_shared_facets] = kfacet;
            ++num_shared_facets;
            break;
          }
        }
      }

      uint num_rvals = sendbuf_facets[dest].size();
      MPI_Sendrecv(&sendbuf_rval[0], num_shared_facets, MPI_UNSIGNED, src, 2,
                   &recvbuf_rval[0], num_rvals, MPI_UNSIGNED, dest, 2,
                   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

      // Add facets with a match on the adjacent rank 'dest'
      for (int k = 0; k < recv_count; ++k)
      {
        shared_facets_.insert(sendbuf_facets[dest][recvbuf_rval[k]]);
      }
    }

    // Cleanup
    delete[] recvbuf_rval;
    delete[] sendbuf_rval;
    delete[] recvbuf_facets;
    delete[] facet_vertices;
    delete[] sendbuf_facets;
    delete[] sendbuf_vertices;

#endif
  }

  message("GlobalFacetMap: number of shared facets %u", shared_facets_.size());
}
//-----------------------------------------------------------------------------
bool GlobalFacetMap::is_global(Facet& facet)
{
  return ((facet.num_entities(tdim_) == 1) &&
          (shared_facets_.count(facet.index()) == 0));
}
//-----------------------------------------------------------------------------
bool GlobalFacetMap::is_shared(Facet& facet)
{
  return ((facet.num_entities(tdim_) == 1) &&
          (shared_facets_.count(facet.index()) > 0));
}
//-----------------------------------------------------------------------------
void GlobalFacetMap::disp() const
{
  section("GlobalFacetMap");
  message("Number of shared facets : %u", shared_facets_.size());
  end();
  skip();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
