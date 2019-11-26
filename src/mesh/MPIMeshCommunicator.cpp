// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/MPIMeshCommunicator.h>

#include <dolfin/common/timing.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(MeshValues<uint, Vertex>& dist)
{

  if (!dist.mesh().is_distributed())
  {
    return;
  }

#if HAVE_MPI

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
      sendbuf_vcoords[owner].append(v->x(), v->x() + gdim);
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
    mesh.swap( new_mesh );
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
  uint * recvbuf_v = local_vgindex.ptr() + local_vgindex_size;
  // Resize vertex coordinates array to fit new cells
  uint const coords_size = coords.size();
  coords.resize(coords_size + recvmax_x);
  real * recvbuf_x = coords.ptr() + coords_size;
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
  mesh.topology().distdata()[0].swap( distdata );
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

#endif /* HAVE_MPI */

}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute( MeshValues< uint, Cell > & dist,
                                      MeshData *                 D )
{

	if ( !dist.mesh().is_distributed() )
	{
		return;
	}

#if HAVE_MPI

	message( 1, "MPIMeshCommunicator : distribute cells" );
	tic();

	Mesh &              mesh    = dist.mesh();
	uint const          pe_rank = mesh.topology().comm_rank();
	uint const          pe_size = mesh.topology().comm_size();
	uint const          tdim    = mesh.topology().dim();
	uint const          gdim    = mesh.geometry().dim();
	DistributedData distdata( mesh.topology().comm() );

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

  // FIXME remove me
  {
    MeshValues< uint, Vertex > v1( mesh, 0 );

    for ( VertexIterator v( mesh ); !v.end(); ++v )
    {
      v1( v->index() ) = v->global_index();

      // if ( v->is_owned() )
      //   v1( *v ) += 1;

      // if ( v->is_ghost() )
      //   v1( *v ) += 2;

      // if ( v->is_shared() )
      //   v1( *v ) += 4;
    }

    File( "before.pvd" ) << v1;
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
  uint vidx = 0;

	// Collect mesh entities according to distribution
  {
    MeshValues< bool, Vertex > vertex_used( mesh, false );

    for ( CellIterator c( mesh ); !c.end(); ++c )
    {
      uint const owner = dist( *c );

      for ( VertexIterator v( *c ); !v.end(); ++v )
      {
        if ( owner != pe_rank )
        {
          send_cells[owner].push_back( v->global_index() );

          if ( not v->is_ghost() and not vertex_used( *v ) )
          {
            vertex_used( *v ) = true;
            send_vgindex[owner].push_back( v->global_index() );
            send_vcoords[owner].append( v->x(), v->x() + gdim );

            // Transfer vertex functions
            if ( not RV.empty() )
            {
              for ( MeshData::iterator< real, Vertex > it( *D ); it.valid(); ++it )
              {
                RV[owner].push_back( ( *it )( *v ) );
              }
            }
          }
        }
        else
        {
          if ( not v->is_ghost() and not vertex_used( *v ) )
          {
            vertex_used( *v ) = true;
            send_vcoords[owner].append( v->x(), v->x() + gdim );
            distdata.set_map( vidx++, v->global_index() );
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

	// Swap local vertex indices
	// Array< uint > local_vgindex;
	// local_vgindex.swap( send_vgindex[pe_rank] );
  // send_vgindex[pe_rank].resize(0);
	dolfin_assert( send_vgindex[pe_rank].empty() );

	// Swap local coordinates
	Array< real > local_vcoords;
	local_vcoords.swap( send_vcoords[pe_rank] );
	dolfin_assert( send_vcoords[pe_rank].empty() );

	// Swap local cells
	Array< uint > local_cells;
	local_cells.swap( send_cells[pe_rank] );
  dolfin_assert( send_cells[pe_rank].empty() );
  dolfin_assert( local_cells.empty() );

	// Exchange the processed entities
	uint recvmax[2] = {0, 0};
	for ( uint j = 0; j < pe_size; ++j )
	{
		uint sendcnt[2] = { static_cast< uint >( send_cells[j].size() ),
		                    static_cast< uint >( send_vgindex[j].size() ) };
		MPI::check_error( MPI_Reduce(
		  sendcnt, recvmax, 2, MPI_UNSIGNED, MPI_SUM, j, distdata.comm() ) );
	}

	uint recvmax_x = recvmax[1] * gdim;

	// Resize cell vertices array to fit new cells
	local_cells.reserve( local_cells.size() + recvmax[0] );

	// Resize vertex indices
	// local_vgindex.reserve( local_vgindex.size() + recvmax[1] );

	// Resize vertex coordinates
	local_vcoords.reserve( local_vcoords.size() + recvmax_x );

	// Naive MeshValues exchange until I fix the ghost bug in the template class
	uint recvmaxUC = recvmax[0] / mesh.type().num_entities( 0 ) * numUC;
	if ( not UC.empty() )
	{
		UC[pe_rank].reserve( UC[pe_rank].size() + recvmaxUC );
	}

	uint recvmaxRV = recvmax[1] * numRV;
	// real *recvbufRV = nullptr;
	if ( not RV.empty() )
	{
		RV[pe_rank].reserve( RV[pe_rank].size() + recvmaxRV );
	}

	// exchange data
	{
		Array< uint > recv_cells( recvmax[0] );
		Array< uint > recv_vgindex( recvmax[1] );
		Array< real > recv_vcoords( recvmax_x );
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
			                          distdata.comm() );

			// Vertices
			n_rcv_v += MPI::sendrecv( send_vgindex[dst].data(),
			                          send_vgindex[dst].size(),
			                          dst,
			                          recv_vgindex.data() + n_rcv_v,
			                          recv_vgindex.size() - n_rcv_v,
			                          src,
			                          1,
			                          distdata.comm() );

			// Coordinates
			n_rcv_x += MPI::sendrecv( send_vcoords[dst].data(),
			                          send_vcoords[dst].size(),
			                          dst,
			                          recv_vcoords.data() + n_rcv_x,
			                          recv_vcoords.size() - n_rcv_x,
			                          src,
			                          2,
			                          distdata.comm() );

			// Transfer cell functions
			if ( not UC.empty() )
			{
				int n = MPI::sendrecv( UC[dst], dst, recv_uc, src, 3, distdata.comm() );
				dolfin_assert( recvmaxUC >= ( uint ) n );
				UC.append( recv_uc.data(), recv_uc.data() + n );
			}

			// Transfer vertex functions
			if ( not RV.empty() )
			{
				int n = MPI::sendrecv( RV[dst], dst, recv_rv, src, 4, distdata.comm() );
				dolfin_assert( recvmaxRV >= ( uint ) n );
				RV.append( recv_rv.data(), recv_rv.data() + n );
			}
		}

		dolfin_assert( n_rcv_v * gdim == n_rcv_x );

		// process (received) vertices
		for ( int i = 0; i < n_rcv_v; ++i )
		{
			if ( not distdata.has_global( recv_vgindex[i] ) )
			{
				distdata.set_map( vidx++, recv_vgindex[i] );
				local_vcoords.append( recv_vcoords.data() + i * gdim,
				                      recv_vcoords.data() + i * gdim + gdim );
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
					if ( not distdata.has_global( v->global_index() ) )
					{
						local_cells.push_back( vidx ); // cl ?!
						local_vcoords.resize( local_vcoords.size() + gdim, 0. );
						distdata.set_map( vidx++, v->global_index() );
						shared_buffer.push_back( v->global_index() );
					}
					else
					{
					  local_cells.push_back( distdata.get_local( v->global_index() ) ); // cl ?!
					}
				}
			}
		}

		// Add new cells
		for ( int i = 0; i < n_rcv_c; ++i )
		{
			if ( distdata.has_global( recv_cells[i] ) )
			{
				local_cells.push_back( distdata.get_local( recv_cells[i] ) ); // cl ?!
			}
			else
			{
				local_cells.push_back( vidx ); // cl ?!
				local_vcoords.resize( local_vcoords.size() + gdim, 0. );
				distdata.set_map( vidx++, recv_cells[i] );
				shared_buffer.push_back( recv_cells[i] );
			}
		}

		// Exchange ghost vertices
		uint shared_count = 0;
		int  r1           = 0;
		int  r2           = 0;
		MPI::all_reduce< MPI::max >( static_cast< uint >( shared_buffer.size() ),
		                             shared_count,
		                             distdata.comm() );
		Array< uint > shared( shared_count );
		Array< real > recv_buff( 100 * shared_buffer.size() * gdim ); // FIXME
		Array< uint > recv_buff_map( 100 * shared_buffer.size() );    // FIXME
		Array< uint > recv_source;
		for ( uint j = 1; j < pe_size; ++j )
		{
			int src = ( pe_rank - j + pe_size ) % pe_size;
			int dst = ( pe_rank + j ) % pe_size;

			// Send ghost vertices to request coordinates
			int recv_count =
			  MPI::sendrecv( shared_buffer, dst, shared, src, 5, distdata.comm() );

			Array< real > send_buff;
			Array< uint > send_buff_indices;

			for ( int j = 0; j < recv_count; j++ )
			{
				if ( distdata.has_global( shared[j] )
			       and not distdata.is_ghost( distdata.get_local( shared[j] ) ) )
				{
					real * local = local_vcoords.data() + distdata.get_local( shared[j] ) * gdim;
					send_buff.append( local, local + gdim );
					send_buff_indices.push_back( shared[j] );

					// if ( not distdata.is_shared( distdata.get_local( shared[j] ) ) )
					// {
					//   distdata.set_shared_adj( distdata.get_local( shared[j] ), src );
					// }
				}
			}

			r1 += MPI::sendrecv( send_buff.data(),
			                     send_buff.size(),
			                     src,
			                     recv_buff.data() + r1,
			                     recv_buff.size() - r1,
			                     dst,
			                     6,
			                     distdata.comm() );
			int r = MPI::sendrecv( send_buff_indices.data(),
			                       send_buff_indices.size(),
			                       src,
			                       recv_buff_map.data() + r2,
			                       recv_buff_map.size() - r2,
			                       dst,
			                       7,
			                       distdata.comm() );

			recv_source.resize( recv_source.size() + r, dst );
			r2 += r;
		}

		dolfin_assert( r1 == r2 * gdim );

		for ( int i = 0; i < r2; ++i )
		{
      if ( not distdata.has_global( recv_buff_map[i] ) )
        distdata.set_map( vidx++, recv_buff_map[i] );

			uint local = distdata.get_local( recv_buff_map[i] );
			distdata.set_ghost( local, recv_source[i] );
			for ( uint j = 0; j < gdim; ++j )
				local_vcoords[local * gdim + j] = recv_buff[i * gdim + j];
		}
	}

	distdata.remap_shared_adj();

  // Clear mesh using swap with new instance
  mesh = Mesh( mesh.type(), mesh.space(), distdata.comm() );

  dolfin_assert( mesh.topology().connectivity( 0 ) == nullptr );

	// NOTE: This implementation only works for homogeneous topologies
	//       Check cell data size just in case.
	if ( ( local_cells.size() % mesh.type().num_entities( 0 ) ) > 0 )
	{
		error( "MPIMeshCommunicator : inconsistent size of cell buffer '%u' (%u)",
		       local_cells.size(), mesh.type().num_entities( 0 ) );
	}
	uint cindex = local_cells.size() / mesh.type().num_entities( 0 );

	// Finalize distributed data
	distdata.finalize();

	// Update topology
	dolfin_assert( vidx == distdata.local_size() );
	mesh.topology().init( 0, vidx );
	mesh.topology().distdata()[0].swap( distdata );
	mesh.topology().init( tdim, cindex );
	mesh.topology()( tdim, 0 ).set( local_cells );
	mesh.topology().finalize();
	dolfin_assert( vidx == mesh.topology().distdata()[0].local_size() );
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
	dolfin_assert( vidx * gdim == local_vcoords.size() );
	mesh.geometry().assign( local_vcoords );
	mesh.geometry().finalize();

	{
		MeshValues< uint, Vertex > v1( mesh, 0 );

		for ( VertexIterator v( mesh ); !v.end(); ++v )
		{
			v1( v->index() ) = v->global_index();
		}

		File( "after.pvd" ) << v1;
	}

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
			it->swap( M );
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
			it->swap( M );
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
			it->swap( M );
		}
	}

	//
	tocd( 1 );

#else
	MAYBE_UNUSED( D );
#endif /* HAVE_MPI */
}
//-----------------------------------------------------------------------------
template<class E>
void MPIMeshCommunicator::check(Mesh& mesh)
{
  if (!mesh.is_distributed()) return;

  uint const tdim = mesh.topology_dimension();
  uint const edim = entity_dimension<E>(mesh);

  if (edim > tdim)
  {
    error("MPIMeshCommunicator : invalid entity dimension %u", edim);
  }

#if HAVE_MPI

  DistributedData& dist = mesh.distdata()[edim];
  uint const pe_rank = dist.comm_rank();
  uint const pe_size = dist.comm_size();

  message("MPIMeshCommunicator : check distribution for dimension %u", edim);

  // Check shared entities adjacency
  {
    Array<uint> * sbuf = new Array<uint> [pe_size];
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
      std::set<uint> global_indices(sbuf[j].begin(), sbuf[j].end());
      if (global_indices.size() != sbuf[j].size())
      {
        error("Duplicate global indices for entities of dimension %u", edim);
      }
    }
    uint const rbuf_size = rbuf.size();
    rbuf.resize(rbuf_size + recv_max);
    uint * recv_buf = rbuf.ptr() + rbuf_size;
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
    delete[] sbuf;
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
      std::set<uint> global_indices(sbuf[j].begin(), sbuf[j].end());
      if (global_indices.size() != sbuf[j].size())
      {
        error("Duplicate global indices for entities of dimension %u", edim);
      }
    }
    uint const rbuf_size = rbuf.size();
    rbuf.resize(rbuf_size + recv_max);
    uint * recv_buf = rbuf.ptr() + rbuf_size;
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
#endif /* HAVE_MPI */
}
//--- TEMPLATE INSTANTIATIONS -------------------------------------------------
template void MPIMeshCommunicator::check<Vertex>(Mesh& Mesh);
template void MPIMeshCommunicator::check<Edge>  (Mesh& Mesh);
template void MPIMeshCommunicator::check<Face>  (Mesh& Mesh);
template void MPIMeshCommunicator::check<Facet> (Mesh& Mesh);
template void MPIMeshCommunicator::check<Cell>  (Mesh& Mesh);
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::check(Mesh& mesh)
{
  uint const tdim = mesh.topology_dimension();
  check<Vertex>(mesh);
  if (tdim > 1) check<Edge>(mesh);
  if (tdim > 2) check<Face>(mesh);
  if (tdim > 0) check<Cell>(mesh);
}
//-----------------------------------------------------------------------------
} /* namespace dolfin */
