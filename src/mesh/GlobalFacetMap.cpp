// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-02-25
// Last changed: 2009-03-03

#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/GlobalFacetMap.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/main/MPI.h>


#ifdef HAVE_MPI
#include <mpi.h>
#endif


namespace dolfin
{

//-----------------------------------------------------------------------------
GlobalFacetMap::GlobalFacetMap(Mesh& mesh) : _mesh(mesh)
{
  init();
}
//-----------------------------------------------------------------------------
GlobalFacetMap::~GlobalFacetMap()
{
}
//-----------------------------------------------------------------------------
#ifdef HAVE_MPI
void GlobalFacetMap::init()
{
  uint const tdim = _mesh.topology().dim();

  // Generate facet - cell connectivity if not generated
  _mesh.init(tdim - 1, tdim);

  // Iterate over all Facets connected to the shared vertices
  for(MeshSharedIterator sv(_mesh.distdata(), 0); !sv.end(); ++sv)
  {
    Vertex v(_mesh, sv.index());
    for(FacetIterator f(v); !f.end(); ++f)
    {
      if (f->num_entities(tdim) == 1)
      {
        // Mark all facets as local facets
        if(global_facet.count(f->index()) == 0)
        {
          global_facet[f->index()] = false;
        }
      }
    }
  }

  switch(tdim)
  {
  case 1:
    findGlobal1D();
    break;
  case 2:
  case 3:
    findGlobalND();
    break;
  default:
    error("Could not handle local to global map with facet of dim %d", tdim);
    break;
  }
}
//-----------------------------------------------------------------------------
void GlobalFacetMap::findGlobal1D()
{
  uint const tdim = _mesh.topology().dim();
  MeshDistributedData const& mddata = _mesh.distdata();

  _map<uint,bool>::iterator iter;
  _mesh.init(tdim - 1, 0);
  for(iter = global_facet.begin(); iter != global_facet.end(); ++iter)
  {
    Vertex v(_mesh, iter->first);

    // Mark as an exterior facet
    if ( v.num_entities(tdim) == 1 && mddata.is_shared(v.index(), 0) )
    {
      global_facet[v.index()] = true;
    }
  }

}
//-----------------------------------------------------------------------------
void GlobalFacetMap::findGlobalND()
{

  Array<uint> send_buff;
  uint const tdim = _mesh.topology().dim();
  MeshDistributedData const& mddata = _mesh.distdata();

  _map<uint,bool>::iterator iter;
  _map<uint, uint>::iterator  uiter;
  _map<uint, uint> unassigned;

  uint num_un = 0;
  _mesh.init(tdim - 1, 0);
  for(iter = global_facet.begin(); iter != global_facet.end(); ++iter)
  {
    Facet f(_mesh, iter->first);
    uint num_shared = 0;

    for (VertexIterator v(f); !v.end(); ++v)
    {
      if( mddata.is_shared(v->index(), 0) )
      {
        ++num_shared;
      }
    }

    if ( f.num_entities(tdim) == 1 && num_shared < f.num_entities(0) )
    {
      global_facet[f.index()] = true;
    }
    else
    {
      unassigned[num_un++] = f.index();
      for (VertexIterator v(f); !v.end(); ++v)
      {
        send_buff.push_back(mddata.get_global(*v));
      }
      // Make an optimistic guess
      global_facet[f.index()] = true;
    }
  }

  int sh_count = send_buff.size();
  //FIXME: Cannot work as it is with heterogeneous mesh since the data packing
  //       is not constant, maybe use the maximum
  uint const num_facet_vertices = _mesh.type().num_vertices(tdim - 1);
  int res_size = sh_count / num_facet_vertices;
  int recv_size = 0;
  int recv_count = 0;

  MPI_Allreduce(&sh_count, &recv_size, 1, MPI_INT,MPI_MAX, MPI::DOLFIN_COMM);
  uint *recv_buff = new uint[ recv_size ];
  uint *res_buff = new uint[ res_size ];

  Array<uint> shared_buff;

  MPI_Status status;
  uint pe_size = MPI::numProcesses();
  uint rank = MPI::processNumber();
  uint src, dest;
  uint num_own = 0;

  for(uint j=1; j<pe_size; ++j)
  {
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, dest, 1,
                 recv_buff, recv_size, MPI_UNSIGNED, src, 1,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);

    for(int i = 0; i < recv_count; i += num_facet_vertices)
    {
      uint vi = 0;
      for(uint k = 0; k < num_facet_vertices; ++k)
      {
        if(! mddata.has_global(recv_buff[i+k], 0))
        {
          ++vi;
        }
      }
      if(vi)
      {
        shared_buff.push_back(0);
        continue;
      }

      Vertex vertex(_mesh, mddata.get_vertex_local(recv_buff[i]));

      for(FacetIterator f(vertex); !f.end(); ++f)
      {
        num_own = 0;

        // Only consider facets connected to one cell
        if ( f->num_entities(tdim) != 1)
        {
          continue;
        }

        for(VertexIterator v(*f); !v.end(); ++v)
        {
          for(uint k = 0; k < tdim; ++k)
          {
            if(recv_buff[i+k] == mddata.get_global(*v))
            {
              ++num_own;
            }
          }
        }
        if(num_own == num_facet_vertices)
        {
          break;
        }
      }
      if(num_own == num_facet_vertices)
      {
        shared_buff.push_back(1);
      }
      else
      {
        shared_buff.push_back(0);
      }

    }

    MPI_Sendrecv(&shared_buff[0], shared_buff.size(), MPI_UNSIGNED, src, 2,
                 res_buff, res_size, MPI_UNSIGNED, dest, 2,
                 MPI::DOLFIN_COMM, &status);

    for(int i = 0; i<res_size; ++i)
    {
      if(res_buff[i] == 1)
      {
        global_facet[unassigned[i]] = false;
      }
    }
    shared_buff.clear();
  }

  delete[] recv_buff;
  delete[] res_buff;
}
//-----------------------------------------------------------------------------
#else
void GlobalFacetMap::init()
{
}
//-----------------------------------------------------------------------------
void GlobalFacetMap::findGlobal1D()
{
}
//-----------------------------------------------------------------------------
void GlobalFacetMap::findGlobalND()
{
}
//-----------------------------------------------------------------------------
#endif
bool GlobalFacetMap::globalFacet(Facet& facet)
{

  const uint index = facet.index();

  // If the facet is in the map, it might be a local facet
  // COMMENT: if the mesh is serial then global_facet is empty then the second
  // test is never evaluated. Conversely the first assertion implies the second
  // as building the global_facet map relies on the shared entities iterator.
  if(global_facet.count(index) > 0 && MPI::numProcesses() > 1)
  {
    return global_facet[index];
  }
  else
  {
    return (facet.num_entities(_mesh.topology().dim()) == 1);
  }
}
//-----------------------------------------------------------------------------

}
