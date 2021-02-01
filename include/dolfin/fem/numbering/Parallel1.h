// Copyright (C) 2013 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRIVATE_NUMBERING_PARALLEL1_H
#define __DOLFIN_PRIVATE_NUMBERING_PARALLEL1_H

#include <dolfin/common/DistributedData.h>
#include <dolfin/fem/DofNumbering.h>

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
 *  @class  Parallel1Numbering
 *
 *  @brief  Implements a parallel dof distribution for vector elements which
 *          ensures that the ownership of dofs matches the ownership of mesh
 *          entities on which they are located.
 *
 */

class Parallel1Numbering : public DofNumbering
{

public:
  ///
  Parallel1Numbering( Mesh & mesh, ufc::dofmap & ufc_dofmap )
    : DofNumbering( mesh, ufc_dofmap )
  {
  }

  ///
  ~Parallel1Numbering() override = default;

  ///
  inline void tabulate_dofs( size_t * dofs,
                             ufc::cell const &,
                             Cell const & cell ) const override
  {
    size_t const ii = ufc_dofmap.num_element_dofs() * cell.index();
    std::copy( &array[ii], &array[ii] + ufc_dofmap.num_element_dofs(), dofs );
  }

  ///
  void build() override
  {
    DofNumbering::init();
    //---
    if ( !mesh.is_distributed() )
    {
      error( "Parallel1Numbering : can only be used for a distributed mesh" );
    }
    //---

#ifdef HAVE_MPI

    /*
     * This is a naive and inefficient renumbering hastly coded, sorry...
     */

    size_t const   tdim = mesh.topology_dimension();
    size_t *       dofs = new size_t[ufc_dofmap.num_element_dofs()];
    _set< size_t > owned;
    _set< size_t > ufc_ghosts;
    _set< size_t > ufc_shared;
    _map< size_t, std::vector< size_t > > dof2index;
    shared_.clear();
    ghosts_.clear();

    size_t pe_size = MPI::size();
    size_t rank    = MPI::rank();

    size_t num_expected_shared = 0;
    size_t num_expected_ghosts = 0;
    for ( size_t i = 0; i < tdim; ++i )
    {
      num_expected_shared +=
        ufc_dofmap.num_entity_dofs( i ) * mesh.topology().num_shared( i );
      num_expected_ghosts +=
        ufc_dofmap.num_entity_dofs( i ) * mesh.topology().num_ghost( i );
    }

    // Cache number of dofs per mesh entity
    size_t * num_entity_dofs    = new size_t[tdim + 1];
    size_t   max_num_entity_dof = 0;
    for ( size_t d = 0; d <= tdim; ++d )
    {
      num_entity_dofs[d] = ufc_dofmap.num_entity_dofs( d );
      max_num_entity_dof = std::max( max_num_entity_dof, num_entity_dofs[d] );
    }
    dolfin_assert( max_num_entity_dof > 0 );
    size_t * entity_local_dofs = new size_t[max_num_entity_dof];
    size_t * entity_dofs       = new size_t[max_num_entity_dof];
    size_t * facet_dofs        = new size_t[ufc_dofmap.num_facet_dofs()];

    size_t       ii = 0;
    CellIterator cell( mesh );
    UFCCell      ufc_cell( *cell );
    for ( ; !cell.end(); ++cell )
    {
      ufc_cell.update( *cell );
      ufc_dofmap.tabulate_dofs( dofs, mesh.num_entities(), ufc_cell.entity_indices );

      // Create mapping from dof to dofmap offset
      for ( size_t i = 0; i < ufc_dofmap.num_element_dofs(); ++i )
      {
        dof2index[dofs[i]].push_back( ii++ );
      }

      // Add dofs restricted to the cell as owned
      size_t const num_celldofs = num_entity_dofs[tdim];
      ufc_dofmap.tabulate_entity_dofs( entity_local_dofs, tdim, 0 );
      for ( size_t dof = 0; dof < num_celldofs; ++dof )
      {
        owned.insert( dofs[entity_local_dofs[dof]] );
      }

      // Decide ownership
      for ( size_t d = 0; d < tdim; ++d )
      {
        size_t const num_dofs = num_entity_dofs[d];
        for ( MeshEntityIterator m( *cell, d ); !m.end(); ++m )
        {
          dolfin_assert( m.pos() == ( size_t ) cell->index( *m ) );
          // Get the dof indices for the entity
          ufc_dofmap.tabulate_entity_dofs( entity_local_dofs, d, m.pos() );
          for ( size_t dof = 0; dof < num_dofs; ++dof )
          {
            entity_dofs[dof] = dofs[entity_local_dofs[dof]];
          }

          //
          if ( m->is_ghost() )
          {
            dolfin_assert( m->is_shared() );
            // Append UFC indices of ghost entities
            ufc_ghosts.insert( entity_dofs, entity_dofs + num_dofs );
          }
          else if ( m->is_shared() )
          {
            // Append owned shared indices: these will not be renumbered
            ufc_shared.insert( entity_dofs, entity_dofs + num_dofs );
            owned.insert( entity_dofs, entity_dofs + num_dofs );
          }
          else
          {
            owned.insert( entity_dofs, entity_dofs + num_dofs );
          }
        }
      }
    }
    delete[] facet_dofs;
    delete[] entity_dofs;
    delete[] entity_local_dofs;
    delete[] num_entity_dofs;
    delete[] dofs;

    dolfin_assert( ufc_ghosts.size() == num_expected_ghosts );
    dolfin_assert( ufc_shared.size()
                   == num_expected_shared - num_expected_ghosts );

    // Renumber dofs
    array_size = mesh.num_cells() * ufc_dofmap.num_element_dofs();
    array      = new size_t[array_size];

    size_t const range       = owned.size();
    size_t       local_index = 0;
    MPI::offset( range, local_index );
    size_t const rank_offset = local_index;
    set_range( rank_offset, range );

    // Compute renumbering for local and owned shared dofs
    size_t                src, dst, max_recv;
    std::vector< size_t > sendbuf;
    for ( _set< size_t >::iterator it = owned.begin(); it != owned.end();
          ++it, ++local_index )
    {
      for ( std::vector< size_t >::iterator di = dof2index[*it].begin();
            di != dof2index[*it].end();
            ++di )
      {
        array[*di] = local_index;
      }
      // Send owned shared dofs: old index then new index
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
    size_t * recvbuf = new size_t[max_recv];
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
        if ( it != ufc_ghosts.end() )
        {
          size_t const newidx = recvbuf[k + 1];
          dolfin_assert( ( newidx < rank_offset )
                         || ( newidx >= rank_offset + range ) );
          ufc_ghosts.erase( it );
          ghosts_.insert( newidx );
          dolfin_assert( shared_.count( newidx ) == 0 );
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

    message(
      1, "Parallel1Numbering: Remaining ghost dofs %u", ufc_ghosts.size() );

    //---

    if ( ghosts_.size() != num_expected_ghosts )
    {
      error( "Mismatch: expected number of shared dofs and actual, %u != %u",
             num_expected_ghosts,
             ghosts_.size() );
    }
    if ( shared_.size() != num_expected_shared )
    {
      error( "Mismatch: expected number of shared dofs and actual, %u != %u",
             num_expected_shared,
             shared_.size() );
    }

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
    return std::string( "Dof numbering for generic parallel vector" );
  }

private:
  ///
  _set< size_t > shared_;
  _set< size_t > ghosts_;
};

}
/* namespace dolfin */

#endif /* __DOLFIN_PRIVATE_NUMBERING_PARALLEL1_H */
