// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/UFCHalo.h>

#include <dolfin/common/AdjacentMapping.h>
#include <dolfin/fem/Coefficient.h>
#include <dolfin/fem/DofMapSet.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/entities/Facet.h>

#include <string>
#include <vector>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFCHalo::UFCHalo( UFC &                                ufc,
                  std::vector< Coefficient * > const & coefficients,
                  DofMapSet const &                    dof_map_set )
  : ufc_( ufc )
  , mesh_( const_cast< Mesh & >( dof_map_set.mesh() ) )
  , coefficients_( coefficients )
  , dof_map_set_( dof_map_set )
  , rank_offsets_()
  , facet_map_()
  , r_packet_size_( 0 )
  , u_packet_size_( 0 )
{
  // Early exit as nothing has to be done.
  if ( !mesh_.is_distributed() )
  {
    return;
  }

  // Clear data structures and define data padding
  clear();
  size_t const      gdim      = mesh_.geometry_dimension();
  size_t const      facet_dim = mesh_.type().facet_dim();
  DistributedData & distdata  = mesh_.distdata()[facet_dim];

  //
  size_t const num_cell_vertices      = mesh_.type().num_entities( 0 );
  size_t const coordinates_data_size  = num_cell_vertices * gdim;
  size_t       coefficients_data_size = 0;

  for ( size_t i = 0; i < ufc_.form.num_coefficients(); ++i )
  {
    coefficients_data_size += ufc_.coefficient_elements[i]->space_dimension();
  }

  size_t dofs_data_size = 0;
  for ( size_t i = 0; i < ufc_.form.rank(); ++i )
  {
    dofs_data_size += ufc_.local_dimensions[i];
  }

  // Separate real and size_t data types to avoid copy of dof indices and casts
  // vertex coordinates + coefficients values
  r_packet_size_ = coordinates_data_size + coefficients_data_size;

  // local facet index + arguments dof indices
  u_packet_size_ = 1 + dofs_data_size;

  // Allocate data structures
  _set< size_t > const & adj               = distdata.get_adj_ranks();
  size_t const           num_shared_facets = distdata.num_shared();

  // Pack by adjacent rank
  size_t offset = 0;
  for ( _set< size_t >::const_iterator it = adj.begin(); it != adj.end(); ++it )
  {
    rank_offsets_.insert( FacetOffsets( *it, offset ) );
    offset += distdata.shared_mapping().to( *it ).size();
  }
  //
  dolfin_assert( offset == num_shared_facets );

  // Cache association of local facet index to offset in halo data
  _map< size_t, size_t > facet_offsets;
  for ( SharedIterator sh( distdata ); sh.valid(); ++sh )
  {
    size_t const ark = *( distdata.get_shared_adj( sh.index() ).begin() );
    size_t       rank_offset = rank_offsets_[ark];

    // Maps local shared ordering to adjacent shared ordering
    std::vector< size_t > const & adjmap = distdata.shared_mapping().to( ark );

    facet_map_[sh.index()] = FacetOffsets( rank_offset + facet_offsets[ark],
                                           rank_offset + adjmap[facet_offsets[ark]] );

    // Increment current facet index for the given adjacent rank
    ++facet_offsets[ark];
  }
  //
  dolfin_assert( facet_map_.size() == num_shared_facets );

  r_data0_.resize( num_shared_facets * r_packet_size_ );
  u_data0_.resize( num_shared_facets * u_packet_size_ );
  r_data1_.resize( num_shared_facets * r_packet_size_ );
  u_data1_.resize( num_shared_facets * u_packet_size_ );

  // Fill data structures
  this->update( coefficients_, dof_map_set_ );
}

//-----------------------------------------------------------------------------

UFCHalo::~UFCHalo()
{
  clear();
}

//-----------------------------------------------------------------------------

void UFCHalo::update( Facet & facet )
{

  FacetMap::const_iterator it = facet_map_.find( facet.index() );
  if ( it == facet_map_.end() )
  {
    if ( facet.is_shared() )
    {
      error( "Shared facet index not found in halo data structure." );
    }
    else
    {
      error( "Trying to fetch halo data for a non-shared facet." );
    }
  }

  // Update real data stored in halo data arrays
  real * r0 = &r_data0_[r_packet_size_ * it->second.first];
  real * r1 = &r_data1_[r_packet_size_ * it->second.second];

  // Update pointers to coordinates
  size_t const gdim = mesh_.geometry_dimension();
  for ( size_t i = 0; i < mesh_.type().num_entities( 0 ); ++i )
  {
    for ( size_t dim = 0; dim < gdim; ++dim, ++r0, ++r1 )
    {
      ufc_.cell0.coordinates[i * Space::MAX_DIMENSION + dim] = *r0;
      ufc_.cell1.coordinates[i * Space::MAX_DIMENSION + dim] = *r1;
    }
  }

  // Update UFC expansion coefficients, needs copy for the moment
  for ( size_t i = 0; i < ufc_.form.num_coefficients(); ++i )
  {
    size_t const spacedim = ufc_.coefficient_elements[i]->space_dimension();
    std::copy( r0, r0 + spacedim, ufc_.macro_w[i] );
    r0 += spacedim;
    std::copy( r1, r1 + spacedim, ufc_.macro_w[i] + spacedim );
    r1 += spacedim;
  }

  // Update size_t data stored in halo data arrays:
  size_t * u0 = &u_data0_[u_packet_size_ * it->second.first];
  size_t * u1 = &u_data1_[u_packet_size_ * it->second.second];

  // Update local facet indices
  ufc_.facet0 = *u0;
  ++u0;
  ufc_.facet1 = *u1;
  ++u1;

  // Update UFC dofs indices for each dimension, needs copy for the moment
  for ( size_t i = 0; i < ufc_.form.rank(); ++i )
  {
    size_t const localdim = ufc_.local_dimensions[i];
    std::copy( u0, u0 + localdim, ufc_.macro_dofs[i] );
    u0 += localdim;
    std::copy( u1, u1 + localdim, ufc_.macro_dofs[i] + localdim );
    u1 += localdim;
  }
}

//-----------------------------------------------------------------------------
void UFCHalo::update( std::vector< Coefficient * > const & coefficients,
                      DofMapSet const &                    dof_map_set )
{
  Mesh & mesh = dof_map_set[0].mesh();

  if ( !mesh.is_distributed() )
  {
    return;
  }

#ifdef HAVE_MPI

  size_t const      tdim      = mesh.topology_dimension();
  size_t const      gdim      = mesh.geometry_dimension();
  size_t const      facet_dim = mesh.type().facet_dim();
  DistributedData & distdata  = mesh.distdata()[facet_dim];

  // Exchange of data for contribution of halo macro elements

  if ( coefficients.size() != ufc_.form.num_coefficients() )
  {
    error( "UFCHalo: invalid number of coefficients passed as arguments:\n"
           "Expected: %d ; Provided: %d",
           ufc_.form.num_coefficients(), coefficients.size() );
  }

  //
  size_t const rank    = MPI::rank();
  size_t const pe_size = MPI::size();

  // Loop over shared facets to collect data following shared iterator ordering
  size_t const num_cell_vertices = mesh_.type().num_entities( 0 );
  for( FacetMap::value_type const & facet_offset : facet_map_ )
  {
    Facet facet( mesh, facet_offset.first );

    // Create cell
    Cell   cell( mesh, facet.entities( tdim )[0] );
    size_t cell_facet = cell.index( facet );
    ufc_.cell.update( cell );

    // Set arrays offset
    real *   r_entry = &r_data0_[facet_offset.second.first * r_packet_size_];
    size_t * u_entry = &u_data0_[facet_offset.second.first * u_packet_size_];

    // Collect data for shared facet contribution
    // TODO: implement proper serialization functions
    for ( size_t i = 0; i < num_cell_vertices; ++i )
    {
      std::copy( &ufc_.cell.coordinates[i * Space::MAX_DIMENSION],
                 &ufc_.cell.coordinates[i * Space::MAX_DIMENSION] + gdim,
                 r_entry );
      r_entry += tdim;
    }

    // Interpolate coefficients on cell
    for ( size_t i = 0; i < ufc_.form.num_coefficients(); ++i )
    {
      coefficients[i]->interpolate( r_entry, ufc_.cell,
                                    *ufc_.coefficient_elements[i],
                                    cell, cell_facet );
      r_entry += ufc_.coefficient_elements[i]->space_dimension();
    }

    // Add local facet index
    *u_entry = cell_facet;
    ++u_entry;

    // Tabulate dofs for each dimension on macro element
    for ( size_t i = 0; i < ufc_.form.rank(); ++i )
    {
      dof_map_set[i].tabulate_dofs( u_entry, ufc_.cell, cell );
      u_entry += ufc_.local_dimensions[i];
    }
  }

  // Exchange data to fill halo data arrays: facet blocks are written directly
  for ( size_t j = 1; j < pe_size; ++j )
  {
    int src = ( rank - j + pe_size ) % pe_size;
    int dst = ( rank + j ) % pe_size;

    size_t const num_send_facets = distdata.shared_mapping().to( dst ).size();
    size_t const num_recv_facets = distdata.shared_mapping().from( src ).size();

    MPI::sendrecv( r_data0_.data() + rank_offsets_[dst] * r_packet_size_,
                   num_send_facets * r_packet_size_, dst,
                   r_data1_.data() + rank_offsets_[src] * r_packet_size_,
                   num_recv_facets * r_packet_size_, src, 1 );

    MPI::sendrecv( u_data0_.data() + rank_offsets_[dst] * u_packet_size_,
                   num_send_facets * u_packet_size_, dst,
                   u_data1_.data() + rank_offsets_[src] * u_packet_size_,
                   num_recv_facets * u_packet_size_, src, 1 );
  }

#else
  MAYBE_UNUSED( coefficients );
#endif
}

//-----------------------------------------------------------------------------
void UFCHalo::disp() const
{
  section( "UFCHalo" );
  prm( "Facet map size", facet_map_.size() );
  prm( "Adjacent ranks", rank_offsets_.size() );
  for ( std::pair< size_t const, size_t > const & offset : rank_offsets_ )
  {
    cout << "\tproc " << offset.first << " : " << offset.second << "\n";
  }
  prm( "Size of real data packet", r_packet_size_ );
  prm( "Size of size_t data packet", u_packet_size_ );
  end();
}

//-----------------------------------------------------------------------------
void UFCHalo::clear()
{
  rank_offsets_.clear();
  facet_map_.clear();
}

} /* namespace dolfin */
