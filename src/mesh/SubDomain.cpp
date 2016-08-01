// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
// Modified by Aurelien Larcher, 2016.
//
// First added:  2007-04-24
// Last changed: 2007-07-21

#include <dolfin/mesh/SubDomain.h>

#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/parameter/parameters.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
SubDomain::SubDomain() :
    abstol_(1.0e-6)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SubDomain::~SubDomain()
{
}
//-----------------------------------------------------------------------------
bool SubDomain::inside(real const * x, bool const on_boundary) const
{
  error("SubDomain : inside() not unimplemented.");
  return false;
}
//-----------------------------------------------------------------------------
bool SubDomain::inside(MeshEntity& entity, bool const on_boundary) const
{
  if (entity.dim() == 0)
  {
    return inside(entity.mesh().geometry().x(entity.index()), on_boundary);
  }
  uint ret = entity.num_entities(0);
  for (VertexIterator v(entity); !v.end(); ++v)
  {
    if (this->inside(v->x(), on_boundary))
    {
      --ret;
    }
  }
  return (ret == 0);
}
//-----------------------------------------------------------------------------
bool SubDomain::overlap(MeshEntity& entity, bool const on_boundary) const
{
  if (entity.dim() == 0)
  {
    return inside(entity.mesh().geometry().x(entity.index()), on_boundary);
  }
  uint ret = false;
  for (VertexIterator v(entity); !v.end(); ++v)
  {
    if (this->inside(v->x(), on_boundary))
    {
      ret = true;
      break;
    }
  }
  return ret;
}
//-----------------------------------------------------------------------------
bool SubDomain::close(real const x, real const xref, real const abstol) const
{
  return (std::fabs(x - xref) < abstol);
}
//-----------------------------------------------------------------------------
bool SubDomain::close(real const x, real const xref) const
{
  return (std::fabs(x - xref) < abstol_);
}
//-----------------------------------------------------------------------------
void SubDomain::mark(MeshFunction<uint>& sub_domains, uint index) const
{
  /*
   message(1, "Computing sub domain markers for sub domain %d.", sub_domain);
   error("WIP");
   */

  Mesh& mesh = sub_domains.mesh();
  uint const tdim = mesh.topology().dim();
  uint const edim = sub_domains.dim();

  // Compute sub domain markers
  bool const is_distributed = mesh.is_distributed();
  for (MeshEntityIterator entity(mesh, edim); !entity.end(); ++entity)
  {
    // Check if entity is on the boundary
    bool on_boundary = false;
    if (edim == tdim - 1)
    {
      on_boundary = (entity->num_entities(tdim) == 1) && !entity->is_shared();
    }
    else if (edim == 0)
    {
      for (FacetIterator fi(*entity); !fi.end(); ++fi)
      {
        if ((fi->num_entities(tdim) == 1) && !fi->is_shared())
        {
          on_boundary = true;
          break;
        }
      }
    }
    //
    if (this->inside(*entity, on_boundary))
    {
      sub_domains(*entity) = index;
    }
  }

#ifdef HAVE_MPI
  if (mesh.is_distributed())
  {
    uint const pe_size = MPI::numProcesses();
    uint const rank = MPI::processNumber();

    Array<uint> * sendbuf = new Array<uint> [pe_size];
    for (GhostIterator it(mesh.distdata()[edim]); !it.end(); ++it)
    {
      if (sub_domains.get(it.index()) == index)
      {
        sendbuf[it.owner()].push_back(it.global_index());
      }
    }

    //
    MPI_Status status;
    uint src;
    uint dst;
    int send_size;
    int recv_size;
    int recv_count;
    for (uint j = 0; j < pe_size; ++j)
    {
      send_size = sendbuf[j].size();
      MPI_Reduce(&send_size, &recv_size, 1, MPI_INT, MPI_SUM, j,
                 MPI::DOLFIN_COMM);
    }
    uint * recvbuf = new uint[recv_size];
    for (uint j = 1; j < pe_size; ++j)
    {
      src = (rank - j + pe_size) % pe_size;
      dst = (rank + j) % pe_size;

      MPI_Sendrecv(&sendbuf[dst][0], sendbuf[dst].size(), MPI_UNSIGNED, dst, 1,
                   &recvbuf[0], recv_size, MPI_UNSIGNED, src, 1,
                   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

      for (int k = 0; k < recv_count; ++k)
      {
        sub_domains.set(mesh.distdata()[edim].get_local(recvbuf[k]), index);
      }
    }
    delete[] recvbuf;
    delete[] sendbuf;
  }
#endif

}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
