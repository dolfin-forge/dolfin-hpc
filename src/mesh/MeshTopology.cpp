// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/MeshTopology.h>

#include <dolfin/log/log.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Connectivity.h>
#include <dolfin/mesh/EntityKey.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/utilities/MeshRenumber.h>

#include <algorithm>
#include <ctime>

namespace dolfin
{

//-----------------------------------------------------------------------------

MeshTopology::MeshTopology( CellType const & type, Comm & comm, bool frozen )
  : Distributed< MeshTopology >( comm )
  , type_( type.clone() )
  , dim_( type.dim() )
  , frozen_( frozen )
  , distdata_( dim_ )
  , timestamp_( 0 )
{
  for ( size_t i = 0; i < MeshTopology::CMAX; ++i )
  {
    for ( size_t j = 0; j < MeshTopology::CMAX; ++j )
    {
      C_[i][j] = nullptr;
    }
  }

  update_token();
}

//-----------------------------------------------------------------------------

MeshTopology::MeshTopology( MeshTopology const & other )
  : Distributed< MeshTopology >( other )
  , type_( cloneptr( other.type_ ) )
  , dim_( other.dim_ )
  , frozen_( other.frozen_ )
  , num_entities_( other.num_entities_ )
  , distdata_( other.distdata_ )
  , timestamp_( other.timestamp_ )
{
  for ( size_t i = 0; i < MeshTopology::CMAX; ++i )
  {
    for ( size_t j = 0; j < MeshTopology::CMAX; ++j )
    {
      C_[i][j] = copyptr( other.C_[i][j] );
    }
  }
}

//-----------------------------------------------------------------------------

MeshTopology::~MeshTopology()
{
  for ( size_t i = 0; i < MeshTopology::CMAX; ++i )
  {
    for ( size_t j = 0; j < MeshTopology::CMAX; ++j )
    {
      delete C_[i][j];
      C_[i][j] = nullptr;
    }
  }

  delete type_;
  type_ = nullptr;
}

//-----------------------------------------------------------------------------

auto swap( MeshTopology & a, MeshTopology & b ) -> void
{
  using std::swap;

  swap( static_cast< Distributed< MeshTopology > & >( a ),
        static_cast< Distributed< MeshTopology > & >( b ) );

  swap( a.type_,         b.type_ );
  swap( a.dim_,          b.dim_ );
  swap( a.frozen_,       b.frozen_ );
  swap( a.num_entities_, b.num_entities_ );
  swap( a.distdata_,     b.distdata_ );
  swap( a.timestamp_,    b.timestamp_ );

  for ( size_t i = 0; i < MeshTopology::CMAX; ++i )
  {
    for ( size_t j = 0; j < MeshTopology::CMAX; ++j )
    {
      swap( a.C_[i][j], b.C_[i][j] );
    }
  }
}

//-----------------------------------------------------------------------------

auto MeshTopology::operator=( const MeshTopology & other ) -> MeshTopology &
{
  MeshTopology tmp( other );
  swap( *this, tmp );

  return *this;
}

//-----------------------------------------------------------------------------

auto MeshTopology::operator==( MeshTopology const & other ) const -> bool
{
  if ( this == &other )
    return true;

  if ( !objptrcmp( type_, other.type_ ) )
    return false;

  for ( size_t i = 0; i <= dim_; ++i )
  {
    if ( this->size( i ) != other.size( i ) )
      return false;
  }

  for ( size_t i = 0; i < MeshTopology::CMAX; ++i )
  {
    for ( size_t j = 0; j < MeshTopology::CMAX; ++j )
    {
      if ( !objptrcmp( C_[i][j], other.C_[i][j] ) )
        return false;
    }
  }

  if ( distdata_ != other.distdata_ )
    return false;

  return true;
}

//-----------------------------------------------------------------------------

auto MeshTopology::operator!=( MeshTopology const & other ) const -> bool
{
  return !( *this == other );
}

//-----------------------------------------------------------------------------

auto MeshTopology::init( size_t dim, size_t nlocal, size_t nglobal ) -> void
{
  // Overflow
  if ( dim_ < dim )
  {
    error( "MeshTopology : initializing entities of dimension %u but topology"
           "dimension is %u",
           dim,
           dim_ );
  }

  dolfin_assert( type_ );
  if ( C_[dim][0] == nullptr )
  {
    dolfin_assert( !C_[dim][0] );
    C_[dim][0] = new Connectivity( nlocal, type_->num_vertices( dim ) );
  }
  else if ( dim_ == 0 )
  {
    // NOTE: point meshes have cell dimension is equal to the vertex dimension
    //       re-initialization is allowed but subject to constraint.
    if ( nlocal != C_[dim][0]->order() )
    {
      error( "MeshTopology : re-initializing point cell with different size" );
    }
  }
  else
  {
    error( "MeshTopology : re-initializing entities of dimension %u", dim );
  }

  // Set size of distributed data
  if ( distributed() )
  {
    if ( nglobal && ( nglobal < nlocal ) )
    {
      error( "MeshTopology : number of global entities lower than number of "
             " local entities %u < %u",
             nlocal,
             nglobal );
    }
    distdata_[dim].set_size( nlocal, nglobal );
  }
  else
  {
    // In serial, require that the number of global entities is not initialized
    // of is equal to the number of local entities
    if ( nglobal && ( nlocal != nglobal ) )
    {
      error( "MeshTopology : invalid number of global entities set in serial" );
    }
  }
}

//-----------------------------------------------------------------------------

auto MeshTopology::finalize() -> void
{
  // Reorder cells according to UFC convention
  reorder();

  // Finalize distributed data
  if ( distributed() )
  {
    for ( size_t d = 0; d <= dim_; ++d )
    {
      if ( this->connectivity( d ) )
      {
        DistributedData & ddata = this->distdata()[d];
        if ( !ddata.is_finalized() )
        {
          ddata.finalize();
        }
        if ( this->size( d ) != ddata.local_size() )
        {
          error( "MeshEditor : vertex size mismatch between topology '%u' and "
                 "distributed data '%u'", this->size( d ), ddata.local_size() );
        }
      }
    }
  }

  // update num_entities_
  this->num_entities_update();

  // Do not renumber automatically !
  // This would cause issues for boundary meshes and some mesh algorithms.
}

//-----------------------------------------------------------------------------

auto MeshTopology::remap( size_t d0, std::vector< size_t > const & mapping )
  -> void
{
  for ( size_t d1 = 0; d1 <= dim_; ++d1 )
  {
    if ( C_[d0][d1] )
    {
      C_[d0][d1]->remap_l( mapping );
      update_token();
    }
    if ( C_[d1][d0] )
    {
      C_[d1][d0]->remap_r( mapping );
      update_token();
    }
  }

  reorder();

  // Remap distributed data
  if ( distributed() )
  {
    distdata_[d0].remap_numbering( mapping );
  }
}

//-----------------------------------------------------------------------------

auto MeshTopology::num_entities_update() const -> void
{
  num_entities_.resize( dim_ + 1 );

  for ( size_t d = 0; d <= dim_; ++d )
  {
    if ( this->connectivity( d ) )
    {
      num_entities_[d] = this->global_size( d );
    }
  }
}

//-----------------------------------------------------------------------------

auto MeshTopology::entities( size_t di ) const -> Connectivity const *
{
  if ( connectivity( di ) )
  {
    return connectivity( di );
  }

  /*
   *  Compute entities to obtain ( tdim, di ) and ( di, 0 ) connectivities
   */

  Connectivity const * const cv = connectivity( dim_ );
  if ( cv != nullptr )
  {
    // Initialize local array of entities
    size_t const                         m = type_->num_entities( di );
    size_t const                         n = type_->num_vertices( di );
    std::vector< std::vector< size_t > > vertex_entities( this->size( 0 ) );
    size_t **                            entities = new size_t *[m];
    for ( size_t e = 0; e < m; ++e )
    {
      entities[e] = new size_t[n]();
    }

    // Collect entities and create cell -> entities connectivities
    EntityKey                  key( n );
    std::vector< EntityKey * > entities_list;
    Connectivity *             ce = new Connectivity( this->size( dim_ ), m );
    for ( size_t c = 0; c < cv->order(); ++c )
    {
      type_->create_entities( entities, di, ( *cv )[c].data() );
      for ( size_t e = 0; e < m; ++e )
      {
        key.set( entities[e], entities_list.size() );
        size_t v0 = entities[e][0];
        dolfin_assert( v0 < this->size( 0 ) );
        if ( vertex_entities[v0].size() == 0 )
        {
          entities_list.push_back( new EntityKey( n, entities[e], key.idx ) );
        }
        else
        {
          // If the first vertex does not contain entity, append new key
          for ( size_t i = 0; i < vertex_entities[v0].size(); ++i )
          {
            if ( *entities_list[vertex_entities[v0][i]] == key )
            {
              key.idx = entities_list[vertex_entities[v0][i]]->idx;
              goto add_entity;
            }
          }
          entities_list.push_back( new EntityKey( n, entities[e], key.idx ) );
        }

      add_entity:
        for ( size_t v = 0; v < n; ++v )
        {
          vertex_entities[entities[e][v]].push_back( key.idx );
        }
        ( *ce )[c][e] = key.idx;
      }
    }
    dolfin_assert( C_[dim_][di] == nullptr );
    std::swap( C_[dim_][di], ce );

    // Create entity -> vertices connectivities from collected entities
    Connectivity * ev = new Connectivity( entities_list.size(), n );
    for ( size_t e = 0; e < entities_list.size(); ++e )
    {
      ev->set( e, entities_list[e]->indices );
    }
    dolfin_assert( C_[di][0] == nullptr );
    std::swap( C_[di][0], ev );

    // Cleanup
    destruct( entities_list );
    for ( size_t e = 0; e < m; ++e )
    {
      delete[] entities[e];
    }
    delete[] entities;
  }
  else
  {
    error( "MeshTopology: computation of entities requires cell vertices" );
  }

  // New entities have been computed: trigger renumbering
  renumber();

  return connectivity( di );
}

//-----------------------------------------------------------------------------

auto MeshTopology::transpose( size_t d1, size_t d0 ) const
  -> Connectivity const *
{
  if ( connectivity( d0, d1 ) )
  {
    return connectivity( d0, d1 );
  }

  Connectivity const * const c10 = connectivity( d1, d0 );

  size_t const          card0 = compute( d0, 0 )->order();
  std::vector< size_t > conn( card0, 0 );

  for ( size_t e = 0; e < c10->order(); ++e )
    for ( size_t i = 0; i < c10->degree( e ); ++i )
      conn[( *c10 )[e][i]]++;

  Connectivity * c01 = new Connectivity( conn );
  std::fill( conn.begin(), conn.end(), 0 );

  for ( size_t e1 = 0; e1 < c10->order(); ++e1 )
  {
    for ( size_t e0 = 0; e0 < c10->degree( e1 ); ++e0 )
    {
      size_t const ei          = ( *c10 )[e1][e0];
      ( *c01 )[ei][conn[ei]++] = e1;
    }
  }
  dolfin_assert( c01->order() == card0 );
  std::swap( C_[d0][d1], c01 );

  return connectivity( d0, d1 );
}
//-----------------------------------------------------------------------------
auto MeshTopology::intersection( size_t d0, size_t di, size_t d1 ) const
  -> Connectivity const *
{
  if ( connectivity( d0, d1 ) )
  {
    return connectivity( d0, d1 );
  }

  /*
   *  Compute connectivities using the intersection d0 - di - d1.
   */

  Connectivity const * const c0v = compute( d0, di );
  Connectivity const * const c1v = compute( d1, di );
  Connectivity const * const cv1 = compute( di, d1 );

  if ( d0 != d1 )
  {
    size_t const                         o0v = c0v->order();
    std::vector< std::vector< size_t > > conn( o0v );
    for ( size_t e0 = 0; e0 < o0v; ++e0 )
    {
      _ordered_set< size_t >        entities;
      std::vector< size_t > const & c0 = ( *c0v )[e0];
      size_t const                  d0 = c0v->degree( e0 );
      for ( size_t i = 0; i < d0; ++i )
      {
        std::vector< size_t > const & c1 = ( *cv1 )[c0[i]];
        size_t const                  d1 = cv1->degree( c0[i] );
        for ( size_t j = 0; j < d1; ++j )
        {
          if ( contains( c0, ( *c1v )[c1[j]] ) )
          {
            entities.insert( c1[j] );
          }
        }
      }
      conn[e0].assign( entities.begin(), entities.end() );
    }
    Connectivity * c01 = new Connectivity( conn );
    dolfin_assert( c01->order() == o0v );
    std::swap( C_[d0][d1], c01 );
  }
  else
  {
    size_t const                         o0v = c0v->order();
    std::vector< std::vector< size_t > > conn( o0v );
    for ( size_t e0 = 0; e0 < o0v; ++e0 )
    {
      _ordered_set< size_t >        entities;
      std::vector< size_t > const & c0 = ( *c0v )[e0];
      size_t const                  d0 = c0v->degree( e0 );
      for ( size_t i = 0; i < d0; ++i )
      {
        std::vector< size_t > const & c1 = ( *cv1 )[c0[i]];
        size_t const                  d1 = cv1->degree( c0[i] );
        for ( size_t j = 0; j < d1; ++j )
        {
          if ( e0 != c1[j] )
          {
            entities.insert( c1[j] );
          }
        }
      }
      conn[e0].assign( entities.begin(), entities.end() );
    }
    Connectivity * c01 = new Connectivity( conn );
    dolfin_assert( c01->order() == o0v );
    std::swap( C_[d0][d1], c01 );
  }

  return connectivity( d0, d1 );
}

//-----------------------------------------------------------------------------

auto MeshTopology::compute( size_t d0, size_t d1 ) const -> Connectivity const *
{
  if ( connectivity( d0, d1 ) )
  {
    return connectivity( d0, d1 );
  }

  if ( connectivity( d1, d0 ) )
  {
    /*
     *  Compute connectivities from transpose if possible.
     *  Requires: ( d1, d0 )
     */
    transpose( d1, d0 );
  }
  else
  {
    size_t const dp = std::max( d0, d1 );
    size_t const dm = std::min( d0, d1 );
    if ( d0 != d1 && ( dp == dim_ || dm == 0 ) )
    {
      /*
       *  Compute entities connectivity ( tdim, di ) and ( di, 0 ).
       *  Requires: ( tdim, 0 )
       */
      entities( dm ? dm : dp );
      compute( d0, d1 );
    }
    else if ( dp )
    {
      /*
       *  Compute connectivities ( d0, d1 ).
       *  Requires: ( d0, 0 ), ( d1, 0 ), ( 0, d1 )
       */
      intersection( dp, 0, dm );
      compute( d0, d1 );
    }
    else
    {
      error( "MeshTopology: vertices were not initialized" );
    }
  }

  /*
   * If the created connectivity needs ordering then trigger reordering
   */
  if ( type_->connectivity_needs_ordering( d0, d1 ) )
  {
    reorder();
  }

  if ( connectivity( d0, d1 ) == nullptr )
  {
    error( "MeshTopology: connectivity (%u, %u) not computed.", d0, d1 );
  }

  // update the number of mesh entities
  this->num_entities_update();

  return connectivity( d0, d1 );
}

//-----------------------------------------------------------------------------

auto MeshTopology::reorder() const -> void
{
  // NOTE: ensure that boundary meshes and empty meshes are not reordered
  if ( not frozen_ and this->connectivity( dim_ ) )
  {
    MeshTopology & topology  = const_cast< MeshTopology & >( *this );
    size_t const   num_cells = this->size( dim_ );
    for ( size_t i = 0; i < num_cells; ++i )
    {
      type_->order_entities( topology, i );
    }
  }
}

//-----------------------------------------------------------------------------

auto MeshTopology::renumber() const -> void
{
  if ( distributed() )
  {
    MeshTopology & topology = const_cast< MeshTopology & >( *this );
    if ( not MeshRenumber::renumber( topology ) )
    {
      warning( "MeshTopology: triggered mesh renumbering for nothing" );
    }

    // update num_entities_
    topology.num_entities_update();
  }
}

//-----------------------------------------------------------------------------

auto MeshTopology::disp() const -> void
{
  section( "MeshTopology" );
  //---
  cout << "Dimension   : " << dim_ << "\n";
  cout << "Distributed : " << distributed() << "\n";
  skip();
  begin( "Number of entities:" );
  for ( size_t d = 1; d <= dim_; ++d )
  {
    if ( C_[d][0] )
    {
      cout << d << ": " << C_[d][0]->order() << "\n";
    }
    else
    {
      cout << d << ": "
           << "x\n";
    }
  }
  end();
  skip();
  begin( "Connectivity:" );
  cout << " ";
  for ( size_t d1 = 0; d1 <= dim_; ++d1 )
  {
    cout << " " << d1;
  }
  cout << "\n";
  for ( size_t d0 = 0; d0 <= dim_; ++d0 )
  {
    cout << d0;
    for ( size_t d1 = 0; d1 <= dim_; ++d1 )
    {
      if ( C_[d0][d1] )
      {
        cout << " x";
      }
      else
      {
        cout << " -";
      }
    }
    cout << "\n";
  }
  cout << "\n";
  end();
  //---
  end();
  skip();
}

//-----------------------------------------------------------------------------

auto MeshTopology::update_token() -> void
{
  timestamp_ = std::time( nullptr );
}

//-----------------------------------------------------------------------------

auto MeshTopology::token() const -> int
{
  return timestamp_ ^ size( 0 ) ^ size( dim_ );
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
