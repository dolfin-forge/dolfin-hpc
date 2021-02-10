// Copyright (C) 2008 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRIVATE_NUMBERING_PARALLEL0_H
#define __DOLFIN_PRIVATE_NUMBERING_PARALLEL0_H

#include <dolfin/fem/DofNumbering.h>

#include <dolfin/common/timing.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>

#include <string>

namespace dolfin
{

/**
 *  @class  Parallel0Numbering
 *
 *  @brief  Implements an optimal randomized parallel dof distribution for
 *          scalar elements which does not ensure that all the dofs are owned
 *          by the owner of the entity on which the node is located.
 *
 */

class Parallel0Numbering : public DofNumbering
{

public:
  ///
  Parallel0Numbering( Mesh & mesh, ufc::dofmap & ufc_dofmap )
    : DofNumbering( mesh, ufc_dofmap )
  {
  }

  ///
  ~Parallel0Numbering() override = default;

  ///
  inline void tabulate_dofs( size_t * dofs,
                             ufc::cell const &,
                             Cell const & cell ) const override
  {
    size_t const ii = ufc_dofmap.num_element_dofs() * cell.index();
    dolfin_assert( array != NULL );
    std::copy( &array[ii], &array[ii] + ufc_dofmap.num_element_dofs(), dofs );
  }

  ///
  void build() override
  {
    DofNumbering::init();
    //---
    if ( !mesh.is_distributed() )
    {
      error( "Parallel0Numbering : can only be used for a distributed mesh" );
    }
    //---

#ifdef HAVE_MPI

    message( 1, "Parallel0Numbering : build" );
    tic();

    _map< size_t, size_t >                dof_vote;
    _map< size_t, std::vector< size_t > > dof2index;
    _set< size_t >                        owned;
    _set< size_t >                        ufc_ghosts;
    _set< size_t >                        ufc_shared;
    shared_.clear();
    ghosts_.clear();

    Cell    c0( mesh, 0 );
    UFCCell ufc_cell( c0 );

    size_t const pe_size = MPI::size();
    size_t const rank    = MPI::rank();

    size_t const tdim          = mesh.topology_dimension();
    size_t const local_dim     = ufc_dofmap.num_element_dofs();
    size_t *     dofs          = new size_t[local_dim];
    size_t const nb_facet_dofs = ufc_dofmap.num_facet_dofs();
    size_t *     facet_dofs    = new size_t[nb_facet_dofs];

    // Initialize random number generator differently on each process
    // FIXME: if the ghosts are randomly generated for a given mesh then
    // two instances for the same mesh and same numbering have different
    // ghosts which leads to several issues:
    // - non matching ghosts with consequence of accessing the wrong dof
    // - different local size which leads to crash
    // We have to remove the randomness until a better solution is found.
    // srand((size_t)time(0) + MPI::rank());
    std::srand( MPI::seed() );

    // List shared dofs and assign vote
    // std::vector<size_t> * sendbuf = new std::vector<size_t>[pe_size];
    std::vector< size_t > sendbuf;
    BoundaryMesh &        interior_boundary = mesh.interior_boundary();
    for ( CellIterator bcell( interior_boundary ); !bcell.end(); ++bcell )
    {
      Facet f( mesh, interior_boundary.facet_index( *bcell ) );
      Cell  c( mesh, f.entities( tdim )[0] );

      ufc_cell.update( c );
      ufc_dofmap.tabulate_dofs( dofs, mesh.topology().num_entities(),
                                ufc_cell.entity_indices );
      size_t local_facet = c.index( f );
      ufc_dofmap.tabulate_facet_dofs( facet_dofs, local_facet );

      for ( size_t i = 0; i < nb_facet_dofs; ++i )
      {
        // Assign an ownership vote for each "shared" dof
        size_t dofidx = dofs[facet_dofs[i]];
        if ( ufc_shared.find( dofidx ) == ufc_shared.end() )
        {
          ufc_shared.insert( dofidx );
          size_t const vote = ( size_t )( std::rand() + ( real ) rank );
          dof_vote[dofidx]  = vote;
          sendbuf.push_back( dofidx );
          sendbuf.push_back( vote );
        }
      }
    }

    // Decide ownership of "shared" dofs
    size_t src;
    size_t dst;
    size_t max_recv = sendbuf.size();
    MPI::all_reduce_in_place< MPI::max >( max_recv );
    size_t * recvbuf = new size_t[max_recv];
    for ( size_t j = 1; j < pe_size; ++j )
    {
      src = ( rank - j + pe_size ) % pe_size;
      dst = ( rank + j ) % pe_size;

      int recv_count = MPI::sendrecv(
        &sendbuf[0], sendbuf.size(), dst, recvbuf, max_recv, src, 1 );

      for ( int i = 0; i < recv_count; i += 2 )
      {
        if ( ufc_shared.find( recvbuf[i] ) != ufc_shared.end() )
        {
          // Move dofs with higher ownership votes from shared to forbidden
          if ( recvbuf[i + 1] < dof_vote[recvbuf[i]]
               || ( recvbuf[i + 1] == dof_vote[recvbuf[i]] && ( src < rank ) ) )
          {
            ufc_ghosts.insert( recvbuf[i] );
            ufc_shared.erase( recvbuf[i] );
          }
        }
      }
    }
    sendbuf.clear();

    // Mark all non forbidden dofs as owned by the processes
    for ( CellIterator cell( mesh ); !cell.end(); ++cell )
    {
      ufc_cell.update( *cell );
      ufc_dofmap.tabulate_dofs( dofs, mesh.topology().num_entities(),
                                ufc_cell.entity_indices );
      for ( size_t i = 0; i < local_dim; ++i )
      {
        if ( ufc_ghosts.find( dofs[i] ) == ufc_ghosts.end() )
        {
          // Mark dof as owned
          owned.insert( dofs[i] );
        }

        // Create mapping from dof to dofmap offset
        dof2index[dofs[i]].push_back( cell->index() * local_dim + i );
      }
    }

    // Cleanup
    delete[] recvbuf;
    delete[] facet_dofs;
    delete[] dofs;

    // Renumber dofs
    array_size = mesh.num_cells() * ufc_dofmap.num_element_dofs();
    array      = new size_t[array_size];

    size_t const range       = owned.size();
    size_t       local_index = 0;
    MPI::offset( range, local_index );
    size_t const rank_offset = local_index;
    set_range( rank_offset, range );

    // Compute renumbering for local and owned shared dofs
    for ( _set< size_t >::iterator it = owned.begin(); it != owned.end();
          ++it, ++local_index )
    {
      for ( std::vector< size_t >::iterator di = dof2index[*it].begin();
            di != dof2index[*it].end();
            ++di )
      {
        array[*di] = local_index;
      }
      if ( ufc_shared.count( *it ) > 0 )
      {
        sendbuf.push_back( *it );
        sendbuf.push_back( local_index );
        shared_.insert( local_index );
      }
    }
    ufc_shared.clear();
    max_recv = sendbuf.size();
    MPI::all_reduce_in_place< MPI::max >( max_recv );
    recvbuf = new size_t[max_recv];
    _set< size_t > new_ghosts;
    for ( size_t j = 1; j < pe_size; ++j )
    {
      src = ( rank - j + pe_size ) % pe_size;
      dst = ( rank + j ) % pe_size;

      int recv_count = MPI::sendrecv(
        &sendbuf[0], sendbuf.size(), dst, recvbuf, max_recv, src, 1 );

      for ( int k = 0; k < recv_count; k += 2 )
      {
        // Assign new dof number for ghost dofs
        _set< size_t >::iterator it = ufc_ghosts.find( recvbuf[k] );
        if ( ufc_ghosts.count( recvbuf[k] ) > 0 )
        {
          size_t const newidx = recvbuf[k + 1];
          dolfin_assert( ( newidx < rank_offset )
                         || ( newidx >= rank_offset + range ) );
          ufc_ghosts.erase( it );
          // Create ghost set with new numbering
          ghosts_.insert( newidx );
          // Append ghost dofs to owned-shared dofs
          shared_.insert( newidx );
          for ( std::vector< size_t >::iterator di =
                  dof2index[recvbuf[k]].begin();
                di != dof2index[recvbuf[k]].end();
                ++di )
          {
            array[*di] = newidx;
          }
        }
      }
    }
    delete[] recvbuf;

    tocd();

#endif /* HAVE_MPI */
  }

  ///
  inline bool is_shared( size_t index ) const override
  {
    return ( shared_.count( index ) > 0 );
  }

  ///
  inline bool is_ghost( size_t index ) const override
  {
    return ( ghosts_.count( index ) > 0 );
  }

  ///
  inline std::string description() const override
  {
    return std::string( "Dof numbering for generic parallel scalar" );
  }

private:
  ///
  _set< size_t > shared_;
  _set< size_t > ghosts_;
};

}
/* namespace dolfin */

#endif /* __DOLFIN_PRIVATE_NUMBERING_PARALLEL0_H */
