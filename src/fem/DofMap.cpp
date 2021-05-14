// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/DofMap.h>

#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>

#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------

namespace helper
{

inline auto init_nedofs( ufc::dofmap const & ufc_dofmap )
  -> std::vector< size_t >
{
  // +1 doesnt hurt here, but it is not the cleanest
  std::vector< size_t > nedofs( ufc_dofmap.topological_dimension() + 1 );

  for ( size_t i = 0; i < nedofs.size(); ++i )
    nedofs[i] = ufc_dofmap.num_entity_dofs( i );

  return nedofs;
}

//-----------------------------------------------------------------------------

inline auto init_necdofs( ufc::dofmap const & ufc_dofmap )
  -> std::vector< size_t >
{
  // +1 doesnt hurt here, but it is not the cleanest
  std::vector< size_t > necdofs( ufc_dofmap.topological_dimension() + 1 );

  for ( size_t i = 0; i < necdofs.size(); ++i )
    necdofs[i] = ufc_dofmap.num_entity_closure_dofs( i );

  return necdofs;
}

//-----------------------------------------------------------------------------

inline auto init_flattened( ufc::dofmap const & ufc_dofmap )
  -> std::vector< ufc::dofmap const * >
{
  std::vector< ufc::dofmap const * > flt;
  DofMap::flatten( &ufc_dofmap, flt );
  return flt;
}

//-----------------------------------------------------------------------------

inline auto init_sddims( ufc::dofmap const & ufc_dofmap )
  -> std::vector< size_t >
{
  std::vector< size_t > sddims;

  // Information for mixed elements
  size_t const nb_sub = ufc_dofmap.num_sub_dofmaps();
  if ( nb_sub > 0 )
  {
    // Set local dimensions
    for ( size_t i = 0; i < nb_sub; ++i )
    {
      ufc::dofmap * subdm = ufc_dofmap.create_sub_dofmap( i );
      sddims.push_back( subdm->num_element_dofs() );
      delete subdm;
    }
  }
  else
  {
    sddims = { ufc_dofmap.num_element_dofs() };
  }

  return sddims;
}

//-----------------------------------------------------------------------------

inline auto init_sdoffs( ufc::dofmap const & ufc_dofmap )
  -> std::vector< size_t >
{
  std::vector< size_t > sdoffs;

  // Information for mixed elements
  size_t const nb_sub = ufc_dofmap.num_sub_dofmaps();
  if ( nb_sub > 0 )
  {
    // Set offsets and local dimensions
    size_t off = 0;
    for ( size_t i = 0; i < nb_sub; ++i )
    {
      ufc::dofmap * subdm = ufc_dofmap.create_sub_dofmap( i );
      sdoffs.push_back( off );
      off += subdm->num_element_dofs();
      delete subdm;
    }
  }
  else
  {
    sdoffs = { 0 };
  }

  return sdoffs;
}

//-----------------------------------------------------------------------------

} // namespace helper

//-----------------------------------------------------------------------------

DofMap::DofMap( Mesh & mesh, ufc::form const & form, size_t const i )
  : MeshDependent( mesh )
  , ufc_dofmap_( form.create_dofmap( i ) )
  , offset_( 0 )
  , numbering_( DofNumbering::create( mesh, *ufc_dofmap_ ) )
  , hash_( make_hash( mesh, *ufc_dofmap_ ) )
  , sub_dofmaps_dims_( helper::init_sddims( *ufc_dofmap_ ) )
  , sub_dofmaps_offs_( helper::init_sdoffs( *ufc_dofmap_ ) )
  , flattened_( helper::init_flattened( *ufc_dofmap_ ) )
  , periodic_dofmap_( nullptr )
  , signature( ufc_dofmap_->signature() )
  , tdim( ufc_dofmap_->topological_dimension() )
  , global_dim( ufc_dofmap_->global_dimension( mesh.topology().num_entities() ) )
  , num_sub_dofmaps( ufc_dofmap_->num_sub_dofmaps() )
  , num_global_support_dofs( ufc_dofmap_->num_global_support_dofs() )
  , num_element_support_dofs( ufc_dofmap_->num_element_support_dofs() )
  , num_element_dofs( ufc_dofmap_->num_element_dofs() )
  , num_facet_dofs( ufc_dofmap_->num_facet_dofs() )
  , num_entity_dofs( helper::init_nedofs( *ufc_dofmap_ ) )
  , num_entity_closure_dofs( helper::init_nedofs( *ufc_dofmap_ ) )
{
  // Build the DOLFIN dofmap
  message( 1, "DofMap: init dofmap for signature:\n %s", signature.c_str() );
  numbering_->build();
  message( 1, "DofMap: offset = %u; size = %u",
           numbering_->offset(), numbering_->size() );
}

//-----------------------------------------------------------------------------

DofMap::DofMap( Mesh & mesh, ufc::dofmap const & dofmap )
  : MeshDependent( mesh )
  , ufc_dofmap_( dofmap.create() )
  , offset_( 0 )
  , numbering_( DofNumbering::create( mesh, *ufc_dofmap_ ) )
  , hash_( make_hash( mesh, *ufc_dofmap_ ) )
  , sub_dofmaps_dims_( helper::init_sddims( *ufc_dofmap_ ) )
  , sub_dofmaps_offs_( helper::init_sdoffs( *ufc_dofmap_ ) )
  , flattened_( helper::init_flattened( *ufc_dofmap_ ) )
  , periodic_dofmap_( nullptr )
  , signature( ufc_dofmap_->signature() )
  , tdim( ufc_dofmap_->topological_dimension() )
  , global_dim( ufc_dofmap_->global_dimension( mesh.topology().num_entities() ) )
  , num_sub_dofmaps( ufc_dofmap_->num_sub_dofmaps() )
  , num_global_support_dofs( ufc_dofmap_->num_global_support_dofs() )
  , num_element_support_dofs( ufc_dofmap_->num_element_support_dofs() )
  , num_element_dofs( ufc_dofmap_->num_element_dofs() )
  , num_facet_dofs( ufc_dofmap_->num_facet_dofs() )
  , num_entity_dofs( helper::init_nedofs( *ufc_dofmap_ ) )
  , num_entity_closure_dofs( helper::init_nedofs( *ufc_dofmap_ ) )
{
  // Build the DOLFIN dofmap
  message( 1, "DofMap: init dofmap for signature:\n %s", signature.c_str() );
  numbering_->build();
  message( 1, "DofMap: offset = %u; size = %u",
           numbering_->offset(), numbering_->size() );
}

//-----------------------------------------------------------------------------

DofMap::DofMap( DofMap const & dofmap, size_t i )
  : MeshDependent( dofmap.mesh() )
  , ufc_dofmap_( dofmap.ufc().create_sub_dofmap( i ) )
  , offset_( 0 )
  , numbering_( DofNumbering::create( this->mesh(), *ufc_dofmap_ ) )
  , hash_( make_hash( this->mesh(), *ufc_dofmap_ ) )
  , sub_dofmaps_dims_( helper::init_sddims( *ufc_dofmap_ ) )
  , sub_dofmaps_offs_( helper::init_sdoffs( *ufc_dofmap_ ) )
  , flattened_( helper::init_flattened( *ufc_dofmap_ ) )
  , periodic_dofmap_( nullptr )
  , signature( ufc_dofmap_->signature() )
  , tdim( ufc_dofmap_->topological_dimension() )
  , global_dim( ufc_dofmap_->global_dimension( mesh().topology().num_entities() ) )
  , num_sub_dofmaps( ufc_dofmap_->num_sub_dofmaps() )
  , num_global_support_dofs( ufc_dofmap_->num_global_support_dofs() )
  , num_element_support_dofs( ufc_dofmap_->num_element_support_dofs() )
  , num_element_dofs( ufc_dofmap_->num_element_dofs() )
  , num_facet_dofs( ufc_dofmap_->num_facet_dofs() )
  , num_entity_dofs( helper::init_nedofs( *ufc_dofmap_ ) )
  , num_entity_closure_dofs( helper::init_nedofs( *ufc_dofmap_ ) )
{
  message(
    1, "DofMap: Extracted dofmap for subspace: %s", ufc_dofmap_->signature() );

  // Build the DOLFIN dofmap
  message( 1, "DofMap: init dofmap for signature:\n %s", signature.c_str() );
  numbering_->build();
  message( 1, "DofMap: offset = %u; size = %u",
           numbering_->offset(), numbering_->size() );
}

//-----------------------------------------------------------------------------

DofMap::DofMap( DofMap const &                dofmap,
                std::vector< size_t > const & subsystem,
                size_t &                      offset )
  : MeshDependent( dofmap.mesh() )
  , ufc_dofmap_( dofmap.create_sub_dofmap( subsystem, offset_ ) )
  , offset_( 0 )
  , numbering_( DofNumbering::create( this->mesh(), *ufc_dofmap_ ) )
  , hash_( make_hash( this->mesh(), *ufc_dofmap_ ) )
  , sub_dofmaps_dims_( helper::init_sddims( *ufc_dofmap_ ) )
  , sub_dofmaps_offs_( helper::init_sdoffs( *ufc_dofmap_ ) )
  , flattened_( helper::init_flattened( *ufc_dofmap_ ) )
  , periodic_dofmap_( nullptr )
  , signature( ufc_dofmap_->signature() )
  , tdim( ufc_dofmap_->topological_dimension() )
  , global_dim( ufc_dofmap_->global_dimension( mesh().topology().num_entities() ) )
  , num_sub_dofmaps( ufc_dofmap_->num_sub_dofmaps() )
  , num_global_support_dofs( ufc_dofmap_->num_global_support_dofs() )
  , num_element_support_dofs( ufc_dofmap_->num_element_support_dofs() )
  , num_element_dofs( ufc_dofmap_->num_element_dofs() )
  , num_facet_dofs( ufc_dofmap_->num_facet_dofs() )
  , num_entity_dofs( helper::init_nedofs( *ufc_dofmap_ ) )
  , num_entity_closure_dofs( helper::init_nedofs( *ufc_dofmap_ ) )
{
  // Check that dof map has not be re-ordered
  offset = offset_;
  message( 1, "DofMap: Extracted dofmap for sub system: %s", signature.c_str() );
  message( 1, "DofMap: Offset for sub system: %d", offset );

  // Reset offset
  offset_ = 0;

  // Build the DOLFIN dofmap
  message( 1, "DofMap: init dofmap for signature:\n %s", signature.c_str() );
  numbering_->build();
  message( 1, "DofMap: offset = %u; size = %u",
           numbering_->offset(), numbering_->size() );
}

//-----------------------------------------------------------------------------

DofMap::DofMap( DofMap const & copy )
  : MeshDependent( copy.mesh() )
  , ufc_dofmap_( copy.ufc_dofmap_->create() )
  , offset_( copy.offset_ )
  , numbering_( DofNumbering::create( this->mesh(), *ufc_dofmap_ ) )
  , hash_( copy.hash_ )
  , sub_dofmaps_dims_( copy.sub_dofmaps_dims_ )
  , sub_dofmaps_offs_( copy.sub_dofmaps_offs_ )
  , flattened_( helper::init_flattened( *ufc_dofmap_ ) )
  , periodic_dofmap_( nullptr )
  , signature( copy.signature )
  , tdim( copy.tdim )
  , global_dim( copy.global_dim )
  , num_sub_dofmaps( copy.num_sub_dofmaps )
  , num_global_support_dofs( copy.num_global_support_dofs )
  , num_element_support_dofs( copy.num_element_support_dofs )
  , num_element_dofs( copy.num_element_dofs )
  , num_facet_dofs( copy.num_facet_dofs )
  , num_entity_dofs( copy.num_entity_dofs )
  , num_entity_closure_dofs( copy.num_entity_closure_dofs )
{
}

//-----------------------------------------------------------------------------

DofMap::~DofMap()
{
  if ( ufc_dofmap_ != nullptr )
    delete ufc_dofmap_;

  if ( numbering_ != nullptr )
    delete numbering_;

  destruct( flattened_ );

  if ( periodic_dofmap_ != nullptr )
    delete periodic_dofmap_;
}

//-----------------------------------------------------------------------------

auto DofMap::periodic_mapping( FiniteElementSpace const & space ) const
  -> PeriodicDofsMapping const &
{
  dolfin_assert( space.dofmap() == *this );

  if ( periodic_dofmap_ == nullptr )
  {
    periodic_dofmap_ = new PeriodicDofsMapping( space );
  }
  return *periodic_dofmap_;
}

//-----------------------------------------------------------------------------

auto DofMap::create_sub_dofmap( std::vector< size_t > const & sub_system ) const
  -> ufc::dofmap *
{
  // Reset offset
  size_t local_offset = 0;

  // Recursively extract sub dof map
  ufc::dofmap * sub_dofmap =
    DofMap::create_sub_dofmap( *ufc_dofmap_, sub_system, local_offset );
  message( 1, "DofMap: Extracted ufc dof map for sub system: %s",
           sub_dofmap->signature() );
  message( 1, "DofMap: Local offset for sub system: %d", local_offset );

  return sub_dofmap;
}

//-----------------------------------------------------------------------------

auto DofMap::create_sub_dofmap( std::vector< size_t > const & sub_system,
                                size_t & local_offset ) const -> ufc::dofmap *
{
  // Reset offset
  local_offset = 0;

  // Recursively extract sub dof map
  ufc::dofmap * sub_dofmap =
    DofMap::create_sub_dofmap( *ufc_dofmap_, sub_system, local_offset );
  message( 1, "DofMap: Extracted ufc dof map for sub system: %s",
           sub_dofmap->signature() );
  message( 1, "DofMap: Local offset for sub system: %d", local_offset );

  return sub_dofmap;
}

//-----------------------------------------------------------------------------

auto DofMap::create_sub_dofmap( ufc::dofmap const &           dofmap,
                                std::vector< size_t > const & sub_system,
                                size_t & local_offset ) -> ufc::dofmap *
{
  // Check that a sub system has been specified
  if ( sub_system.size() == 0 )
  {
    // error("Unable to extract sub system (no sub system specified).");
    return dofmap.create();
  }

  // Check if there are any sub systems
  if ( dofmap.num_sub_dofmaps() == 0 )
  {
    error( "Unable to extract sub system (there are no sub systems)." );
  }

  // Check the number of available sub systems
  if ( sub_system[0] >= dofmap.num_sub_dofmaps() )
  {
    error( "Unable to extract sub system %d (only %d sub systems defined).",
           sub_system[0], dofmap.num_sub_dofmaps() );
  }

  // Add to offset if necessary
  for ( size_t i = 0; i < sub_system[0]; i++ )
  {
    ufc::dofmap * ufc_sub_dofmap = dofmap.create_sub_dofmap( i );
    local_offset += ufc_sub_dofmap->num_element_dofs();
    delete ufc_sub_dofmap;
  }

  // Create sub system
  ufc::dofmap * sub_dofmap = dofmap.create_sub_dofmap( sub_system[0] );

  // Return sub system if sub sub system should not be extracted
  if ( sub_system.size() == 1 )
    return sub_dofmap;

  // Otherwise, recursively extract the sub sub system
  std::vector< size_t > sub_sub_system;
  for ( size_t i = 1; i < sub_system.size(); i++ )
  {
    sub_sub_system.push_back( sub_system[i] );
  }
  ufc::dofmap * sub_sub_dofmap =
    DofMap::create_sub_dofmap( *sub_dofmap, sub_sub_system, local_offset );
  delete sub_dofmap;

  return sub_sub_dofmap;
}

//-----------------------------------------------------------------------------

void DofMap::flatten( ufc::dofmap const *                  dofmap,
                      std::vector< ufc::dofmap const * > & stack )
{
  // Single root element or max level is set to zero, return immediately
  if ( dofmap->num_sub_dofmaps() == 0 )
  {
    stack.push_back( dofmap->create() );
    return;
  }
  // Go one level down
  for ( size_t s = 0; s < dofmap->num_sub_dofmaps(); ++s )
  {
    ufc::dofmap const * sub = dofmap->create_sub_dofmap( s );
    if ( sub->num_sub_dofmaps() == 0 )
    {
      // Leaf dofmap
      stack.push_back( sub );
    }
    else
    {
      // Branch
      DofMap::flatten( sub, stack );
    }
  }
}

//-----------------------------------------------------------------------------

auto DofMap::can_vectorize( std::vector< ufc::dofmap const * > flattened )
  -> bool
{
  for ( size_t s = 1; s < flattened.size(); ++s )
  {
    if ( std::strcmp( flattened[0]->signature(), flattened[s]->signature() )
         != 0 )
    {
      return false;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------

void DofMap::disp() const
{
  section( "DofMap" );
  begin( "ufc::dofmap info" );
  prm( "Signature",              signature );
  prm( "Topological dimension",  tdim );
  prm( "Global dimension",       global_dim );
  prm( "# global support dofs",  num_global_support_dofs );
  prm( "# element support dofs", num_element_support_dofs );
  prm( "# element dofs",         num_element_dofs );
  prm( "# facet dofs",           num_facet_dofs );
  prm( "# entity dofs",          to_string( num_entity_dofs ) );
  prm( "# entity closure dofs",  to_string( num_entity_closure_dofs ) );
  end();
  end();
}

//-----------------------------------------------------------------------------

auto DofMap::check( bool throw_error ) -> bool
{
  bool ret = true;

  message( "Check dofs distribution" );
  message( "signature   : %s", signature.c_str() );
  bool distributed_by_entities_ = true;
  message( "by entities : %d", distributed_by_entities_ );
  // FIXME Check if a dof is owned twice, send to everyone to make sure
  BoundaryMesh &        boundary = mesh().interior_boundary();
  Mesh &                mesh     = this->mesh();
  MeshDistributedData & distdata = mesh.distdata();
  size_t const          tdim     = mesh.topology_dimension();

  if ( ufc_dofmap_->num_facet_dofs() == 0 )
  {
    return true;
  }

  //
  _ordered_set< size_t >                                      shared_owned;
  typedef _ordered_map< size_t, std::pair< size_t, size_t > > EntitiesDofMap;
  EntitiesDofMap shared_owned_entities;
  Cell           c0( mesh, 0 );
  UFCCell        ufc_cell( c0 );
  size_t *       cell_dofs       = new size_t[ufc_dofmap_->num_element_dofs()];
  size_t *       num_entity_dofs = new size_t[tdim];
  for ( size_t d = 0; d < tdim; ++d )
  {
    mesh.init( tdim, d );
    mesh.init( mesh.type().facet_dim(), d );
    num_entity_dofs[d] = ufc_dofmap_->num_entity_dofs( d );
  }
  size_t const num_facet_dofs  = ufc_dofmap_->num_facet_dofs();
  size_t *     loc_entity_dofs = new size_t[ufc_dofmap_->num_facet_dofs()];
  for ( CellIterator bcell( boundary ); !bcell.end(); ++bcell )
  {
    Facet f( boundary.mesh(), boundary.facet_index( *bcell ) );
    dolfin_assert( f.is_shared() );

    // Tabulate cell dofs
    Cell cell( mesh, f.entities( tdim )[0] );
    ufc_cell.update( cell );
    this->tabulate_dofs( cell_dofs, ufc_cell );

    // Check that all shared facet dofs are shared
    size_t const local_facet = cell.index( f );
    ufc_dofmap_->tabulate_facet_dofs( loc_entity_dofs, local_facet );
    _ordered_set< size_t > shared_facet_dofs;
    _ordered_set< size_t > local_facet_dofs;
    for ( size_t dof = 0; dof < num_facet_dofs; ++dof )
    {
      size_t gdof = cell_dofs[loc_entity_dofs[dof]];
      if ( !is_shared( gdof ) )
      {
        error( "Dof with local_index %u on shared facet is not shared",
               loc_entity_dofs[dof] );
      }
      shared_facet_dofs.insert( gdof );
      local_facet_dofs.insert( loc_entity_dofs[dof] );
    }

    // Get boundary facet dofs
    size_t const facet_dim = mesh.type().facet_dim();
    ufc_dofmap_->tabulate_entity_dofs(
      loc_entity_dofs, facet_dim, local_facet );
    for ( size_t dof = 0; dof < num_entity_dofs[facet_dim]; ++dof )
    {
      dolfin_assert( is_shared( cell_dofs[loc_entity_dofs[dof]] ) );
      if ( !is_ghost( cell_dofs[loc_entity_dofs[dof]] ) )
      {
        if ( distributed_by_entities_ )
        {
          if ( distdata[f.dim()].is_ghost( f.index() ) )
          {
            error( "Non-ghosted dof's facet is ghosted" );
          }
          shared_owned_entities[cell_dofs[loc_entity_dofs[dof]]] =
            std::pair< size_t, size_t >(
              distdata[f.dim()].get_global( f.index() ), f.dim() );
        }
        shared_owned.insert( cell_dofs[loc_entity_dofs[dof]] );
      }
    }

    // Check lower dimensional entities of boundary facet
    for ( size_t d = 0; d < ( facet_dim ); ++d )
    {
      for ( MeshEntityIterator m( f, d ); !m.end(); ++m )
      {
        dolfin_assert( m->is_shared() );
        // Get the dof indices for the entity
        size_t     local_index = cell.index( *m );
        MeshEntity e( mesh, d, cell.entities( d )[local_index] );
        if ( !distdata[e.dim()].is_shared( e.index() ) )
        {
          error( "Entity of dim %u in shared facet is not shared", d );
        }
        ufc_dofmap_->tabulate_entity_dofs( loc_entity_dofs, d, local_index );

        for ( size_t dof = 0; dof < num_entity_dofs[d]; ++dof )
        {
          if ( local_facet_dofs.count( loc_entity_dofs[dof] ) == 0 )
          {
            error( "Dof index is not a facet dof:\n"
                   "Entity of dim %u with local index %u",
                   d, cell.index( *m ) );
          }
          size_t gdof = cell_dofs[loc_entity_dofs[dof]];
          if ( shared_facet_dofs.count( gdof ) == 0 )
          {
            error( "Dof index is not a shared facet dof:\n"
                   "Entity of dim %u with local index %u",
                   d, cell.index( *m ) );
          }
          dolfin_assert( is_shared( gdof ) );
          if ( !is_ghost( gdof ) )
          {
            if ( distributed_by_entities_ )
            {
              if ( distdata[m->dim()].is_ghost( m->index() ) )
              {
                error( "Non-ghosted dof's entity of dim %u is ghosted", d );
              }
              shared_owned_entities[gdof] = std::pair< size_t, size_t >(
                distdata[m->dim()].get_global( m->index() ), m->dim() );
            }
            shared_owned.insert( gdof );
          }
        }
      }
    }
  }
  delete[] num_entity_dofs;
  delete[] loc_entity_dofs;
  delete[] cell_dofs;

#if DOLFIN_HAVE_MPI

  // Simple version
  size_t rank    = dolfin::MPI::rank();
  size_t pe_size = dolfin::MPI::size();

  std::vector< size_t > sendbuf;
  if ( distributed_by_entities_ )
  {
    for ( EntitiesDofMap::const_iterator it = shared_owned_entities.begin();
          it != shared_owned_entities.end(); ++it )
    {
      sendbuf.push_back( it->second.first );
      sendbuf.push_back( it->second.second );
      sendbuf.push_back( it->first );
    }
  }
  else
  {
#ifdef __SUNPRO_CC
    for ( _ordered_set< size_t >::iterator it = shared_owned.begin();
          it != shared_owned.end(); ++it )
    {
      sendbuf.push_back( *it );
    }
#else
    sendbuf.insert( sendbuf.end(), shared_owned.begin(), shared_owned.end() );
#endif
  }

  int maxrev_count = sendbuf.size();
  MPI::all_reduce_in_place< MPI::max >( maxrev_count );
  dolfin_assert( maxrev_count > 0 );
  std::vector< size_t > recvbuf( maxrev_count );

  _set< size_t > owned_more;
  for ( size_t p = 1; p < pe_size; ++p )
  {
    int src  = ( rank - p + pe_size ) % pe_size;
    int dest = ( rank + p ) % pe_size;

    //
    int recv_count = MPI::sendrecv( sendbuf, dest, recvbuf, src, 0 );
    for ( size_t i = 0; i < size_t( recv_count ); )
    {
      if ( distributed_by_entities_ )
      {
        size_t const idx = recvbuf[i++];
        size_t const dim = recvbuf[i++];
        if ( distdata[dim].has_global( idx ) )
        {
          size_t loc_id = distdata[dim].get_local( idx );
          if ( !distdata[dim].is_ghost( loc_id ) )
          {
            if ( throw_error )
            {
              error( "Entity on rank %u is owned by rank %u", rank, src );
            }
            else
            {
              warning( "Entity on rank %u is owned by rank %u", rank, src );
            }
          }
        }
      }
      size_t const dof = recvbuf[i++];
      if ( is_shared( dof ) && !is_ghost( dof ) )
      {
        owned_more.insert( dof );
        if ( throw_error )
        {
          error( "Degree of freedom on rank %u is also owned by rank %u",
                 rank, src );
        }
        else
        {
          warning( "Degree of freedom on rank %u is also owned by rank %u",
                   rank, src );
        }
      }
    }
  }
  ret &= owned_more.empty();
#else
  MAYBE_UNUSED( throw_error );
#endif
  return ret;
}

//-----------------------------------------------------------------------------

} // namespace dolfin
