// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/SubDomain.h>

#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/EdgeIterator.h>
#include <dolfin/mesh/entities/iterators/FaceIterator.h>
#include <dolfin/mesh/entities/iterators/FacetIterator.h>
#include <dolfin/parameter/parameters.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
SubDomain::SubDomain()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
template < class Entity >
void SubDomain::mark( MeshValues< size_t, Entity > & sub_domains,
                      size_t                         index ) const
{
  message( 1, "SubDomain: Computing markers for sub domain %d.", index );

  Mesh & mesh = sub_domains.mesh();

  // Compute sub domain markers
  for ( typename Entity::iterator e( mesh ); !e.end(); ++e )
  {
    if ( this->enclosed( *e, e->on_boundary() ) )
    {
      sub_domains( *e ) = index;
    }
  }

#ifdef DOLFIN_HAVE_MPI
  if ( mesh.is_distributed() )
  {
    size_t const      pe_size  = MPI::size();
    size_t const      pe_rank  = MPI::rank();
    DistributedData & distdata = mesh.distdata()[sub_domains.dim()];

    std::vector< std::vector< size_t > > sendbuf( pe_size );

    // Update entities to adjacent ranks.
    // The previous implementation updates only ghost to the owner, which
    // assumes that data will only be used by the owner, maybe not.
    for ( SharedIterator it( distdata ); it.valid(); ++it )
    {
      if ( sub_domains( it.index() ) == index )
      {
        it.adj_enqueue( sendbuf, it.global_index() );
      }
    }

    //
    int recv_size = 0;
    for ( size_t j = 0; j < pe_size; ++j )
    {
      int send_size = sendbuf[j].size();
      MPI::reduce< MPI::sum >( &send_size, &recv_size, 1, j, distdata.comm() );
    }
    std::vector< size_t > recvbuf( recv_size );
    for ( size_t j = 1; j < pe_size; ++j )
    {
      int src = ( pe_rank - j + pe_size ) % pe_size;
      int dst = ( pe_rank + j ) % pe_size;

      int recv_count =
        MPI::sendrecv( sendbuf[dst], dst, recvbuf, src, 1, distdata.comm() );

      for ( int k = 0; k < recv_count; ++k )
      {
        sub_domains( distdata.get_local( recvbuf[k] ) ) = index;
      }
    }
  }
#endif
}

//--- TEMPLATE INSTANTIATIONS -------------------------------------------------

// template bool SubDomain::enclosed(Vertex& entity, bool on_boundary) const;
template bool SubDomain::enclosed( Edge & entity, bool on_boundary ) const;
template bool SubDomain::enclosed( Face & entity, bool on_boundary ) const;
template bool SubDomain::enclosed( Facet & entity, bool on_boundary ) const;
template bool SubDomain::enclosed( Cell & entity, bool on_boundary ) const;

// template bool SubDomain::overlap(Vertex& entity, bool on_boundary) const;
template bool SubDomain::overlap( Edge & entity, bool on_boundary ) const;
template bool SubDomain::overlap( Face & entity, bool on_boundary ) const;
template bool SubDomain::overlap( Facet & entity, bool on_boundary ) const;
template bool SubDomain::overlap( Cell & entity, bool on_boundary ) const;

template void SubDomain::mark( MeshValues< size_t, Vertex > & sub_domains,
                               size_t                         index ) const;
template void SubDomain::mark( MeshValues< size_t, Edge > & sub_domains,
                               size_t                       index ) const;
template void SubDomain::mark( MeshValues< size_t, Face > & sub_domains,
                               size_t                       index ) const;
template void SubDomain::mark( MeshValues< size_t, Facet > & sub_domains,
                               size_t                        index ) const;
template void SubDomain::mark( MeshValues< size_t, Cell > & sub_domains,
                               size_t                       index ) const;

//-----------------------------------------------------------------------------

} /* namespace dolfin */
