// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/UFCHalo.h>

#include <dolfin/common/AdjacentMapping.h>
#include <dolfin/fem/Coefficient.h>
#include <dolfin/fem/Form.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/entities/Facet.h>

#include <string>
#include <vector>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFCHalo::UFCHalo( Form & form )
  : facet0( 0 )
  , facet1( 0 )
  , cell0( Cell( form.mesh(), 0 ) )
  , cell1( Cell( form.mesh(), 0 ) )
  , form_( form )
  , rank_offsets_()
  , facet_map_()
  , r_packet_size_( 0 )
  , u_packet_size_( 0 )
{
  // Early exit as nothing has to be done.
  if ( not form_.mesh().is_distributed() )
  {
    return;
  }

  // Clear data structures and define data padding
  clear();
  size_t const      gdim      = form_.mesh().geometry_dimension();
  size_t const      facet_dim = form_.mesh().type().facet_dim();
  DistributedData & distdata  = form_.mesh().distdata()[facet_dim];

  //
  size_t const num_cell_vertices      = form_.mesh().type().num_entities( 0 );
  size_t const coordinates_data_size  = num_cell_vertices * gdim;
  size_t       coefficients_data_size = 0;

  for ( size_t i = 0; i < form_.num_coefficients(); ++i )
  {
    coefficients_data_size += form.elements()[form.rank() + i].space_dim;
  }

  size_t dofs_data_size = 0;
  for ( size_t i = 0; i < form_.rank(); ++i )
  {
    dofs_data_size += form_.dofmaps()[i]->num_element_dofs;
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
  update();
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
  size_t const gdim = form_.mesh().geometry_dimension();
  for ( size_t i = 0; i < form_.mesh().type().num_entities( 0 ); ++i )
  {
    for ( size_t dim = 0; dim < gdim; ++dim, ++r0, ++r1 )
    {
      cell0.coordinates[i * Space::MAX_DIMENSION + dim] = *r0;
      cell1.coordinates[i * Space::MAX_DIMENSION + dim] = *r1;
    }
  }

  // Update UFC expansion coefficients, needs copy for the moment
  for ( size_t i = 0; i < form_.num_coefficients(); ++i )
  {
    size_t const spacedim = form_.elements()[form_.rank() + i].space_dim;
    std::copy( r0, r0 + spacedim, form_.cache().macro_w[i] );
    r0 += spacedim;
    std::copy( r1, r1 + spacedim, form_.cache().macro_w[i] + spacedim );
    r1 += spacedim;
  }

  // Update size_t data stored in halo data arrays:
  size_t * u0 = &u_data0_[u_packet_size_ * it->second.first];
  size_t * u1 = &u_data1_[u_packet_size_ * it->second.second];

  // Update local facet indices
  facet0 = *u0;
  ++u0;
  facet1 = *u1;
  ++u1;

  // Update UFC dofs indices for each dimension, needs copy for the moment
  for ( size_t i = 0; i < form_.rank(); ++i )
  {
    size_t const localdim = form_.dofmaps()[i]->num_element_dofs;
    std::copy( u0, u0 + localdim, form_.cache().macro_dofs[i] );
    u0 += localdim;
    std::copy( u1, u1 + localdim, form_.cache().macro_dofs[i] + localdim );
    u1 += localdim;
  }
}

//-----------------------------------------------------------------------------
void UFCHalo::update()
{
  Mesh & mesh = form_.mesh();

  if ( !mesh.is_distributed() )
  {
    return;
  }

#ifdef HAVE_MPI

  size_t const      tdim      = mesh.topology_dimension();
  size_t const      facet_dim = mesh.type().facet_dim();
  DistributedData & distdata  = mesh.distdata()[facet_dim];

  // Exchange of data for contribution of halo macro elements

  //
  size_t const rank    = MPI::rank();
  size_t const pe_size = MPI::size();

  // Loop over shared facets to collect data following shared iterator ordering
  size_t const num_cell_vertices = mesh.type().num_entities( 0 );
  for( FacetMap::value_type const & facet_offset : facet_map_ )
  {
    Facet facet( mesh, facet_offset.first );

    // Create cell
    Cell   cell( mesh, facet.entities( tdim )[0] );
    size_t cell_facet = cell.index( facet );
    form_.cache().cell.update( cell );

    // Set arrays offset
    real *   r_entry = &r_data0_[facet_offset.second.first * r_packet_size_];
    size_t * u_entry = &u_data0_[facet_offset.second.first * u_packet_size_];

    // Collect data for shared facet contribution
    // TODO: implement proper serialization functions
    size_t const cell_gdim = form_.cache().cell.geometric_dimension;
    for ( size_t i = 0; i < num_cell_vertices; ++i )
    {
      std::copy( &form_.cache().cell.coordinates[i * cell_gdim],
                 &form_.cache().cell.coordinates[i * cell_gdim] + cell_gdim,
                 r_entry );
      r_entry += tdim;
    }

    // Interpolate coefficients on cell
    for ( size_t i = 0; i < form_.num_coefficients(); ++i )
    {
      FiniteElement const & fe = form_.elements()[form_.rank() + i];
      form_.coefficients()[i]->interpolate( r_entry, form_.cache().cell,
                                            fe.ufc(), cell, cell_facet );
      r_entry += fe.space_dim;
    }

    // Add local facet index
    *u_entry = cell_facet;
    ++u_entry;

    // Tabulate dofs for each dimension on macro element
    for ( size_t i = 0; i < form_.rank(); ++i )
    {
      form_.dofmaps()[i]->tabulate_dofs( u_entry, form_.cache().cell, cell );
      u_entry += form_.dofmaps()[i]->num_element_dofs;
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
