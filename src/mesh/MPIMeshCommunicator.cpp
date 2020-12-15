// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/MPIMeshCommunicator.h>

#include <dolfin/common/timing.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/CellIterator.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/VertexIterator.h>

#include <algorithm>

namespace dolfin
{

namespace MPIMeshCommunicator
{

//-----------------------------------------------------------------------------
void distribute(MeshValues<uint, Vertex>& dist)
{

  if (!dist.mesh().is_distributed())
  {
    return;
  }

#if DOLFIN_HAVE_MPI

  message(1, "MPIMeshCommunicator : distribute vertices");
  tic();

  Mesh& mesh = dist.mesh();
  uint const pe_rank = mesh.topology().comm_rank();
  uint const pe_size = mesh.topology().comm_size();
  uint const tdim = mesh.topology().dim();
  uint const gdim = mesh.geometry().dim();

  // Save global number of vertices to check consistency
  dolfin_assert(mesh.topology().connectivity(0));
  dolfin_assert(mesh.topology().distdata()[0].is_finalized());
  uint const num_global_vertices = mesh.topology().global_size(0);
  if (mesh.topology().size(0) != dist.size())
  {
    error("MPIMeshCommunicator : mismatch between number of vertices and size "
          "of the distribution");
  }
  for (uint d = 1; d <= tdim; ++d)
  {
    if (mesh.topology().connectivity(d) && mesh.topology().size(d))
    {
      error("MPIMeshCommunicator : distribution by vertices but entities of "
            "dimension %u exist", d);
    }
  }

  Array<uint> * sendbuf_vgindex = new Array<uint> [pe_size];
  Array<real> * sendbuf_vcoords = new Array<real> [pe_size];

  // Collect mesh entities according to distribution
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    if (v->is_owned())
    {
      uint const owner = dist(*v);
      append(sendbuf_vcoords[owner],v->x(), v->x() + gdim);
      sendbuf_vgindex[owner].push_back(v->global_index());
    }
  }

  // Swap local indices
  Array<uint> local_vgindex; local_vgindex.swap(sendbuf_vgindex[pe_rank]);
  dolfin_assert(sendbuf_vgindex[pe_rank].size() == 0);

  // Swap local coordinates
  Array<real> coords; coords.swap(sendbuf_vcoords[pe_rank]);
  dolfin_assert(sendbuf_vcoords[pe_rank].size() == 0);

  // Create distributed data
  DistributedData distdata(mesh.distdata()[0].comm());

  // Clear mesh using swap with new instance
  {
    Mesh new_mesh(mesh.type(), mesh.space(), distdata.comm());
    swap( mesh, new_mesh );
  }
  dolfin_assert(mesh.topology().connectivity(0) == nullptr);
  dolfin_assert(mesh.topology().connectivity(tdim) == nullptr);

  // Exchange the vertices
  uint recvmax_v;
  for (uint j = 0; j < pe_size; ++j)
  {
    uint s = sendbuf_vgindex[j].size();
    MPI::check_error( MPI_Reduce(&s, &recvmax_v, 1, MPI_UNSIGNED, MPI_SUM, j,
                                 distdata.comm()) );
  }
  uint recvmax_x = recvmax_v * gdim;
  // Resize vertex indices
  uint const local_vgindex_size = local_vgindex.size();
  local_vgindex.resize(local_vgindex_size + recvmax_v);
  uint * recvbuf_v = local_vgindex.data() + local_vgindex_size;
  // Resize vertex coordinates array to fit new cells
  uint const coords_size = coords.size();
  coords.resize(coords_size + recvmax_x);
  real * recvbuf_x = coords.data() + coords_size;
  int recv_cellsount;
  for (uint j = 1; j < pe_size; ++j)
  {
    int src = (pe_rank - j + pe_size) % pe_size;
    int dst = (pe_rank + j) % pe_size;

    // Vertices
    recv_cellsount = MPI::sendrecv( sendbuf_vgindex[dst].data(), sendbuf_vgindex[dst].size(), dst,
                                &recvbuf_v[0], recvmax_v, src,
                                0, distdata.comm() );
    recvbuf_v += recv_cellsount;
    recvmax_v -= recv_cellsount;

    // Coordinates
    recv_cellsount = MPI::sendrecv( sendbuf_vcoords[dst].data(), sendbuf_vcoords[dst].size(), dst,
                                recvbuf_x, recvmax_x, src,
                                1, distdata.comm() );
    recvbuf_x += recv_cellsount;
    recvmax_x -= recv_cellsount;

  }

  // Cleanup and finalize distributed data
  delete[] sendbuf_vcoords;
  delete[] sendbuf_vgindex;

  // Map local vertex indices
  distdata.set_map(local_vgindex);

  // Finalize distributed data
  distdata.finalize();

  // Update topology
  dolfin_assert(local_vgindex.size() == distdata.local_size());
  mesh.topology().init(0 , distdata.local_size());
  swap( mesh.topology().distdata()[0], distdata );
  mesh.topology().finalize();
  dolfin_assert(local_vgindex.size() == mesh.topology().distdata()[0].local_size());
  if(num_global_vertices != mesh.topology().global_size(0))
  {
    error("MPIMeshCommunicator : vertex distribution :\n"
          "invalid global number of vertices %u != %u",
          num_global_vertices, mesh.topology().global_size(0));
  }

  // Update geometry
  dolfin_assert(local_vgindex.size() * gdim == coords.size());
  mesh.geometry().assign(coords);
  mesh.geometry().finalize();

  //
  tocd(1);

#endif /* DOLFIN_HAVE_MPI */

}
//-----------------------------------------------------------------------------
void distribute( MeshValues< uint, Cell > & dist, MeshData * D )
{
  if ( !dist.mesh().is_distributed() )
  {
    return;
  }

#if DOLFIN_HAVE_MPI

  message( 1, "MPIMeshCommunicator : distribute cells" );
  tic();

	Mesh & mesh = dist.mesh();
	DistributedData distdata( mesh.topology().comm() );
	MPI::Communicator & comm = mesh.topology().comm();

	uint const pe_rank = mesh.topology().comm_rank();
	uint const pe_size = mesh.topology().comm_size();
	uint const tdim    = mesh.topology().dim();
	uint const gdim    = mesh.geometry().dim();

	// Save global number of vertices and cells to check consistency
  dolfin_assert( mesh.topology().connectivity( 0 ) );
  dolfin_assert( mesh.topology().distdata()[0].is_finalized() );
  uint const num_global_vertices = mesh.topology().global_size( 0 );

  dolfin_assert( mesh.topology().connectivity( tdim ) );
  dolfin_assert( mesh.topology().distdata()[tdim].is_finalized() );
  uint const num_global_cells = mesh.topology().global_size( tdim );

  if ( mesh.topology().size( tdim ) != dist.size() )
  {
    error( "MPIMeshCommunicator : mismatch between number of cells and size of "
           "the distribution" );
  }

  Array< Array< uint > > send_cells( pe_size );
  Array< Array< uint > > send_vgindex( pe_size );
  Array< Array< real > > send_vcoords( pe_size );

  /// Support the same cell and vertex function as before
  uint const             numRV = D != nullptr ? D->size< real, Vertex >() : 0;
  Array< Array< real > > RV( ( numRV > 0 ) ? pe_size : 0 );
  uint const             numUC =
    D != nullptr ? ( D->size< bool, Cell >() + D->size< uint, Cell >() ) : 0;
  Array< Array< uint > > UC( ( numUC > 0 ) ? pe_size : 0 );

  /// Check that only desired data is present in D
  if ( D != nullptr && ( D->size() != ( numRV + numUC ) ) )
  {
    error( "MPIMeshCommunicator : transferring <real, Vertex>, <bool, Cell>, "
           "and <uint, Cell> only is supported." );
  }

  // vertex index counter
  uint local_vindex = 0;

  // Collect mesh entities according to distribution
  {
    MeshValues< bool, Vertex > vertex_used( mesh, false );

    for ( CellIterator c( mesh ); !c.end(); ++c )
    {
      uint const owner = dist( *c );

      if ( owner != pe_rank )
      {
        for ( VertexIterator v( *c ); !v.end(); ++v )
        {
          send_cells[owner].push_back( v->global_index() );

          if ( not vertex_used( *v ) and v->is_owned() )
          {
            vertex_used( *v ) = true;
            send_vgindex[owner].push_back( v->global_index() );
            for ( uint l = 0; l < gdim; ++l )
              send_vcoords[owner].push_back( v->x()[l] );

            // Transfer vertex functions
            if ( not RV.empty() )
            {
              for ( MeshData::iterator< real, Vertex > it( *D ); it.valid();
                    ++it )
              {
                RV[owner].push_back( ( *it )( *v ) );
              }
            }
          }
        }
      }
      else
      {
        for ( VertexIterator v( *c ); !v.end(); ++v )
        {
          if ( not vertex_used( *v ) and v->is_owned() )
          {
            vertex_used( *v ) = true;
            for ( uint l = 0; l < gdim; ++l )
              send_vcoords[owner].push_back( v->x()[l] );
            distdata.set_map( local_vindex++, v->global_index() );
          }
        }
      }

      // Transfer cell functions
      if ( not UC.empty() )
      {
        for ( MeshData::iterator< bool, Cell > it( *D ); it.valid(); ++it )
        {
          UC[owner].push_back( ( *it )( *c ) );
        }

        for ( MeshData::iterator< uint, Cell > it( *D ); it.valid(); ++it )
        {
          UC[owner].push_back( ( *it )( *c ) );
        }
      }
    }
  }

  // local vertex indices (our own vertex indices should be empty until now)
  dolfin_assert( send_vgindex[pe_rank].empty() );

  // Swap local coordinates
  Array< real > local_vcoords;
  local_vcoords.swap( send_vcoords[pe_rank] );
  dolfin_assert( send_vcoords[pe_rank].empty() );

  // local cells (our own cells should be empty until now)
  dolfin_assert( send_cells[pe_rank].empty() );
  Array< uint > local_cells;

  // Exchange the processed entities
  uint recvmax[2] = {0, 0};
  for ( uint j = 0; j < pe_size; ++j )
  {
    uint sendcnt[2] = { static_cast< uint >( send_cells[j].size() ),
                        static_cast< uint >( send_vgindex[j].size() ) };
    MPI::check_error( MPI_Reduce(
      sendcnt, recvmax, 2, MPI_UNSIGNED, MPI_SUM, j, comm ) );
  }

  // Resize cell vertices array to fit new cells
  local_cells.reserve( local_cells.size() + recvmax[0] );

  // Resize vertex coordinates
  local_vcoords.reserve( local_vcoords.size() + recvmax[1] * gdim );

  // Naive MeshValues exchange until I fix the ghost bug in the template class
  uint recvmaxUC = recvmax[0] / mesh.type().num_entities( 0 ) * numUC;
  if ( not UC.empty() )
  {
    UC[pe_rank].reserve( UC[pe_rank].size() + recvmaxUC );
  }

  uint recvmaxRV = recvmax[1] * numRV;
  if ( not RV.empty() )
  {
    RV[pe_rank].reserve( RV[pe_rank].size() + recvmaxRV );
  }

  // exchange data
  {
    Array< uint > recv_cells( recvmax[0] );
    Array< uint > recv_vgindex( recvmax[1] );
    Array< real > recv_vcoords( recvmax[1] * gdim );
    Array< uint > recv_uc( recvmaxUC );
    Array< real > recv_rv( recvmaxRV );

    int n_rcv_c = 0;
    int n_rcv_v = 0;
    int n_rcv_x = 0;

    for ( uint j = 1; j < pe_size; ++j )
    {
      int src = ( pe_rank - j + pe_size ) % pe_size;
      int dst = ( pe_rank + j ) % pe_size;

      // Cells
      n_rcv_c += MPI::sendrecv( send_cells[dst].data(),
                                send_cells[dst].size(),
                                dst,
                                recv_cells.data() + n_rcv_c,
                                recv_cells.size() - n_rcv_c,
                                src,
                                0,
                                comm );

      // Vertices
      n_rcv_v += MPI::sendrecv( send_vgindex[dst].data(),
                                send_vgindex[dst].size(),
                                dst,
                                recv_vgindex.data() + n_rcv_v,
                                recv_vgindex.size() - n_rcv_v,
                                src,
                                1,
                                comm );

      // Coordinates
      n_rcv_x += MPI::sendrecv( send_vcoords[dst].data(),
                                send_vcoords[dst].size(),
                                dst,
                                recv_vcoords.data() + n_rcv_x,
                                recv_vcoords.size() - n_rcv_x,
                                src,
                                2,
                                comm );

      // Transfer cell functions
      if ( not UC.empty() )
      {
        int n = MPI::sendrecv( UC[dst], dst, recv_uc, src, 3, comm );
        dolfin_assert( recvmaxUC >= ( uint ) n );
        append( UC[pe_rank], recv_uc.data(), recv_uc.data() + n );
      }

      // Transfer vertex functions
      if ( not RV.empty() )
      {
        int n = MPI::sendrecv( RV[dst], dst, recv_rv, src, 4, comm );
        dolfin_assert( recvmaxRV >= ( uint ) n );
        append( RV[pe_rank], recv_rv.data(), recv_rv.data() + n );
      }
    }

    dolfin_assert( n_rcv_v * gdim == static_cast<uint>( n_rcv_x ) );

    // process (received) vertices
    for ( int i = 0; i < n_rcv_v; ++i )
    {
      if ( not distdata.has_global( recv_vgindex[i] ) )
      {
        distdata.set_map( local_vindex++, recv_vgindex[i] );
        for ( uint l = 0; l < gdim; ++l )
          local_vcoords.push_back( recv_vcoords[i * gdim + l] );
      }
    }

    // process (recveived) cells
    Array< uint > shared_buffer;
    // Add old cells
    for ( CellIterator c( mesh ); !c.end(); ++c )
    {
      if ( dist( *c ) == pe_rank )
      {
        for ( VertexIterator v( *c ); !v.end(); ++v )
        {
          if ( distdata.has_global( v->global_index() ) )
          {
            local_cells.push_back( distdata.get_local( v->global_index() ) );
          }
          else
          {
            local_cells.push_back( local_vindex );
            for ( uint l = 0; l < gdim; ++l )
              local_vcoords.push_back( DOLFIN_REAL_MIN );
            distdata.set_map( local_vindex++, v->global_index() );
            shared_buffer.push_back( v->global_index() );
          }
        }
      }
    }

    // Add new cells
    for ( int i = 0; i < n_rcv_c; ++i )
    {
      if ( distdata.has_global( recv_cells[i] ) )
      {
        local_cells.push_back( distdata.get_local( recv_cells[i] ) );
      }
      else
      {
        local_cells.push_back( local_vindex );
        for ( uint l = 0; l < gdim; ++l )
          local_vcoords.push_back( DOLFIN_REAL_MIN );
        distdata.set_map( local_vindex++, recv_cells[i] );
        shared_buffer.push_back( recv_cells[i] );
      }
    }

    // Exchange ghost vertices
    uint shared_count = shared_buffer.size();
    MPI::all_reduce_in_place< MPI::max >( shared_count, comm );
    Array< uint > shared( shared_count );
    Array< real > recv_buff( shared_buffer.size() * gdim );
    Array< uint > recv_buff_map( shared_buffer.size() );

    for ( uint j = 1; j < pe_size; ++j )
    {
      int src = ( pe_rank - j + pe_size ) % pe_size;
      int dst = ( pe_rank + j ) % pe_size;

      // Send ghost vertices to request coordinates
      int n_recv = MPI::sendrecv( shared_buffer, dst, shared, src, 5, comm );

      Array< real > send_buff;
      Array< uint > send_buff_indices;

      for ( int i = 0; i < n_recv; ++i )
      {
        if ( distdata.has_global( shared[i] ) )
        {
          uint local = distdata.get_local( shared[i] );

          if ( not distdata.is_ghost( local ) )
          {
            bool valid = true;
            for ( uint l = 0; l < gdim; ++l )
              if ( local_vcoords[local * gdim + l] == DOLFIN_REAL_MIN )
                valid = false;

            if ( valid )
            {
              for ( uint l = 0; l < gdim; ++l )
                send_buff.push_back( local_vcoords[local * gdim + l] );

              send_buff_indices.push_back( shared[i] );
              distdata.set_shared_adj( local, src );
            }
          }
        }
      }

      int n_rcv_X = MPI::sendrecv( send_buff, src, recv_buff, dst, 6, comm );
      int n_rcv_v = MPI::sendrecv( send_buff_indices, src, recv_buff_map, dst, 7, comm );

      dolfin_assert( static_cast< uint >( n_rcv_X ) == n_rcv_v * gdim );

      for ( int i = 0; i < n_rcv_v; ++i )
      {
        if ( distdata.has_global( recv_buff_map[i] ) )
        {
          uint local = distdata.get_local( recv_buff_map[i] );

          for ( uint l = 0; l < gdim; ++l )
              local_vcoords[local * gdim + l] = recv_buff[i * gdim + l];

          distdata.set_ghost( local, dst );
        }
      }
    }
  }

  distdata.remap_shared_adj();

#if defined( DEBUG )
  distdata.check_ghost();
  distdata.check_shared();
#endif

  // Finalize distributed data
  distdata.finalize();

  // Clear mesh using swap with new instance
  mesh = Mesh( mesh.type(), mesh.space(), comm );

  dolfin_assert( mesh.topology().connectivity( 0 ) == nullptr );

  // NOTE: This implementation only works for homogeneous topologies
  //       Check cell data size just in case.
  if ( ( local_cells.size() % mesh.type().num_entities( 0 ) ) > 0 )
  {
    error( "MPIMeshCommunicator : inconsistent size of cell buffer '%u' (%u)",
           local_cells.size(), mesh.type().num_entities( 0 ) );
  }
  uint cindex = local_cells.size() / mesh.type().num_entities( 0 );

  // Update topology
  dolfin_assert( local_vindex == distdata.local_size() );
  mesh.topology().init( 0, local_vindex );
  swap( mesh.topology().distdata()[0], distdata );
  mesh.topology().init( tdim, cindex );
  mesh.topology()( tdim, 0 ).set( local_cells );
  mesh.topology().finalize();
  dolfin_assert( local_vindex == mesh.topology().distdata()[0].local_size() );
  if ( num_global_vertices != mesh.topology().global_size( 0 ) )
  {
    error( "MPIMeshCommunicator : cell distribution :\n"
           "invalid global number of vertices %u != %u",
           num_global_vertices,
           mesh.topology().global_size( 0 ) );
  }
  if ( num_global_cells != mesh.topology().global_size( tdim ) )
  {
    error( "MPIMeshCommunicator : cell distribution :\n"
           "invalid global number of cells %u != %u",
           num_global_cells,
           mesh.topology().global_size( tdim ) );
  }

  // Update geometry
  dolfin_assert( local_vindex * gdim == local_vcoords.size() );
  mesh.geometry().resize( local_vindex ); // ?!
  mesh.geometry().assign( local_vcoords );
  mesh.geometry().finalize();

  mesh.extent_update();

  // Recreate mesh functions
  if ( not UC.empty() )
  {
    Array< uint > & mUC( UC[pe_rank] );

    uint ii = 0;
    for ( MeshData::iterator< bool, Cell > it( *D ); it.valid(); ++it, ++ii )
    {
      MeshValues< bool, Cell > M( mesh );
      uint const               nUC = mUC.size();
      dolfin_assert( nUC == M.size() * numUC );
      for ( uint j = ii, k = 0; j < nUC; j += numUC, ++k )
      {
        M( k ) = mUC[j];
      }
      swap( *it, M );
    }
    for ( MeshData::iterator< uint, Cell > it( *D ); it.valid(); ++it, ++ii )
    {
      MeshValues< uint, Cell > M( mesh );
      uint const               nUC = mUC.size();
      dolfin_assert( nUC == M.size() * numUC );
      for ( uint j = ii, k = 0; j < nUC; j += numUC, ++k )
      {
        M( k ) = mUC[j];
      }
      swap( *it, M );
    }
  }

  if ( not RV.empty() )
  {
    Array< real > & mRV( RV[pe_rank] );

    uint ii = 0;
    for ( MeshData::iterator< real, Vertex > it( *D ); it.valid(); ++it, ++ii )
    {
      MeshValues< real, Vertex > M( mesh );
      uint const                 nRV = mRV.size();
      dolfin_assert( nRV
                     == ( M.size() - mesh.distdata()[0].num_ghost() ) * numRV );
      for ( uint j = ii, k = 0; j < nRV; j += numRV, ++k )
      {
        M( k ) = mRV[j];
      }
      swap( *it, M );
    }
  }

  //
  tocd( 1 );

#else
  MAYBE_UNUSED( D );
#endif /* DOLFIN_HAVE_MPI */
}
//-----------------------------------------------------------------------------
template<class E>
void check(Mesh& mesh)
{
  if (!mesh.is_distributed()) return;

  uint const tdim = mesh.topology_dimension();
  uint const edim = entity_dimension<E>(mesh);

  if (edim > tdim)
  {
    error("MPIMeshCommunicator : invalid entity dimension %u", edim);
  }

#if DOLFIN_HAVE_MPI

  DistributedData& dist = mesh.distdata()[edim];
  uint const pe_rank = dist.comm_rank();
  uint const pe_size = dist.comm_size();

  message("MPIMeshCommunicator : check distribution for dimension %u", edim);

  // Check shared entities adjacency
  {
    Array< Array<uint> > sbuf( pe_size );
    uint e_count = 0;
    for (typename E::shared e(mesh); e.valid(); ++e, ++e_count)
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
    uint recv_max;
    for (uint j = 0; j < pe_size; ++j)
    {
      uint s = sbuf[j].size();
      MPI::check_error( MPI_Reduce(&s, &recv_max, 1, MPI_UNSIGNED, MPI_SUM, j,
                                   dist.comm()) );
      // Check that no duplicate global entity was added
      _ordered_set<uint> global_indices(sbuf[j].begin(), sbuf[j].end());
      if (global_indices.size() != sbuf[j].size())
      {
        error("Duplicate global indices for entities of dimension %u", edim);
      }
    }
    uint const rbuf_size = rbuf.size();
    rbuf.resize(rbuf_size + recv_max);
    uint * recv_buf = rbuf.data() + rbuf_size;
    int recv_cellsount;
    for (uint j = 1; j < pe_size; ++j)
    {
      int src = (pe_rank - j + pe_size) % pe_size;
      int dst = (pe_rank + j) % pe_size;

      recv_cellsount = MPI::sendrecv( sbuf[dst].data(), sbuf[dst].size(), dst,
                                  &recv_buf[0], recv_max, src,
                                  0, dist.comm() );

      for (int k = 0; k < recv_cellsount; ++k)
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

      recv_buf += recv_cellsount;
      recv_max -= recv_cellsount;
    }
  }

  // Check ghost entities adjacency
  {
    Array<uint> * sbuf = new Array<uint> [pe_size];
    uint e_count = 0;
    for (typename E::ghost e(mesh); e.valid(); ++e, ++e_count)
    {
      sbuf[e.owner()].push_back(e.global_index());
    }
    // Check that ghost entities indices were counted correctly
    dolfin_assert(e_count == mesh.topology().num_ghost(edim));

    // Swap local entities to array
    Array<uint> rbuf; rbuf.swap(sbuf[pe_rank]);
    dolfin_assert(sbuf[pe_rank].size() == 0);

    // Exchange entities
    uint recv_max;
    for (uint j = 0; j < pe_size; ++j)
    {
      uint s = sbuf[j].size();
      MPI::check_error( MPI_Reduce(&s, &recv_max, 1, MPI_UNSIGNED, MPI_SUM, j,
                                   dist.comm()) );
      // Check that no duplicate global entity was added
      _ordered_set<uint> global_indices(sbuf[j].begin(), sbuf[j].end());
      if (global_indices.size() != sbuf[j].size())
      {
        error("Duplicate global indices for entities of dimension %u", edim);
      }
    }
    uint const rbuf_size = rbuf.size();
    rbuf.resize(rbuf_size + recv_max);
    uint * recv_buf = rbuf.data() + rbuf_size;
    int recv_cellsount;
    _set<uint> recv_idx;
    for (uint j = 1; j < pe_size; ++j)
    {
      int src = (pe_rank - j + pe_size) % pe_size;
      int dst = (pe_rank + j) % pe_size;

      recv_cellsount = MPI::sendrecv( sbuf[dst].data(), sbuf[dst].size(), dst,
                                  &recv_buf[0], recv_max, src,
                                  0, dist.comm() );

      for (int k = 0; k < recv_cellsount; ++k)
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

      recv_buf += recv_cellsount;
      recv_max -= recv_cellsount;
    }

    delete[] sbuf;
  }
#endif /* DOLFIN_HAVE_MPI */
}
//--- TEMPLATE INSTANTIATIONS -------------------------------------------------
template void check<Vertex>(Mesh& Mesh);
template void check<Edge>  (Mesh& Mesh);
template void check<Face>  (Mesh& Mesh);
template void check<Facet> (Mesh& Mesh);
template void check<Cell>  (Mesh& Mesh);
//-----------------------------------------------------------------------------
void check(Mesh& mesh)
{
  uint const tdim = mesh.topology_dimension();
  check<Vertex>(mesh);
  if (tdim > 1) check<Edge>(mesh);
  if (tdim > 2) check<Face>(mesh);
  if (tdim > 0) check<Cell>(mesh);
}
//-----------------------------------------------------------------------------

} /* namespace MPIMeshCommunicator */

} /* namespace dolfin */
