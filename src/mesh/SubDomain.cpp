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
template <>
bool SubDomain::enclosed(Vertex& entity, bool on_boundary) const
{
  return inside(entity.x(), on_boundary);
}
//-----------------------------------------------------------------------------
template <class Entity>
bool SubDomain::enclosed(Entity& entity, bool on_boundary) const
{
  for (VertexIterator v(entity); !v.end(); ++v)
  {
    if (!this->inside(v->x(), on_boundary))
    {
      return false;
    }
  }
  return true;
}
//-----------------------------------------------------------------------------
template <>
bool SubDomain::overlap(Vertex& entity, bool on_boundary) const
{
  return inside(entity.x(), on_boundary);
}
//-----------------------------------------------------------------------------
template <class Entity>
bool SubDomain::overlap(Entity& entity, bool on_boundary) const
{
  for (VertexIterator v(entity); !v.end(); ++v)
  {
    if (this->inside(v->x(), on_boundary))
    {
      return true;
    }
  }
  return false;
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
template <class Entity>
void SubDomain::mark(MeshValues<uint, Entity>& sub_domains, uint index) const
{
  message(1, "Computing sub domain markers for sub domain %d.", index);

  Mesh& mesh = sub_domains.mesh();

  // Compute sub domain markers
  for (typename Entity::iterator e(mesh); !e.end(); ++e)
  {
    if (this->enclosed(*e, e->on_boundary())) { sub_domains(*e) = index; }
  }

#ifdef HAVE_MPI
  if (mesh.is_distributed())
  {
    uint const pe_size = MPI::size();
    uint const pe_rank = MPI::rank();
    DistributedData& distdata = mesh.distdata()[sub_domains.dim()];

    Array<uint> * sendbuf = new Array<uint> [pe_size];

    // Update entities to adjacent ranks.
    // The previous implementation updates only ghost to the owner, which
    // assumes that data will only be used by the owner, maybe not.
    for (SharedIterator it(distdata); it.valid(); ++it)
    {
      if (sub_domains(it.index()) == index)
      {
        it.adj_enqueue(sendbuf, it.global_index());
      }
    }

    //
    int send_size;
    int recv_size;
    for (uint j = 0; j < pe_size; ++j)
    {
      send_size = sendbuf[j].size();
      MPI::check_error( MPI_Reduce(&send_size, &recv_size, 1, MPI_INT, MPI_SUM,
                                   j, distdata.comm()) );
    }
    uint * recvbuf = (recv_size ? new uint[recv_size] : NULL);
    for (uint j = 1; j < pe_size; ++j)
    {
      int src = (pe_rank - j + pe_size) % pe_size;
      int dst = (pe_rank + j) % pe_size;

      int recv_count = MPI::sendrecv( &sendbuf[dst][0], sendbuf[dst].size(), dst,
                                      &recvbuf[0], recv_size, src,
                                      1, distdata.comm() );

      for (int k = 0; k < recv_count; ++k)
      {
        sub_domains(distdata.get_local(recvbuf[k])) = index;
      }
    }

    delete[] recvbuf;
    delete[] sendbuf;
  }
#endif

}

//--- TEMPLATE INSTANTIATIONS -------------------------------------------------

// template bool SubDomain::enclosed(Vertex& entity, bool on_boundary) const;
template bool SubDomain::enclosed(Edge  & entity, bool on_boundary) const;
template bool SubDomain::enclosed(Face  & entity, bool on_boundary) const;
template bool SubDomain::enclosed(Facet & entity, bool on_boundary) const;
template bool SubDomain::enclosed(Cell  & entity, bool on_boundary) const;

// template bool SubDomain::overlap(Vertex& entity, bool on_boundary) const;
template bool SubDomain::overlap(Edge  & entity, bool on_boundary) const;
template bool SubDomain::overlap(Face  & entity, bool on_boundary) const;
template bool SubDomain::overlap(Facet & entity, bool on_boundary) const;
template bool SubDomain::overlap(Cell  & entity, bool on_boundary) const;

template void SubDomain::mark(MeshValues<uint, Vertex>& sub_domains, uint index) const;
template void SubDomain::mark(MeshValues<uint, Edge>  & sub_domains, uint index) const;
template void SubDomain::mark(MeshValues<uint, Face>  & sub_domains, uint index) const;
template void SubDomain::mark(MeshValues<uint, Facet> & sub_domains, uint index) const;
template void SubDomain::mark(MeshValues<uint, Cell>  & sub_domains, uint index) const;

//-----------------------------------------------------------------------------

} /* namespace dolfin */
