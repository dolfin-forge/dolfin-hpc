// Copyright (C) 2016-2017 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/Connectivity.h>

#include <dolfin/log/LogStream.h>
#include <dolfin/log/log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

Connectivity::Connectivity( size_t order, size_t degree )
  : order_( order )
  , min_degree_( degree )
  , max_degree_( degree )
  , connections_( order_, std::vector< size_t >( degree, 0 ) )
{
}

//-----------------------------------------------------------------------------

Connectivity::Connectivity( std::vector< size_t > const & valency )
  : order_( valency.size() )
  , min_degree_( ( order_ > 0 ) ? valency[0] : 0 )
  , max_degree_( 0 )
  , connections_( order_ )
{
  for ( size_t e = 0; e < order_; ++e )
  {
    min_degree_     = std::min( min_degree_, valency[e] );
    max_degree_     = std::max( max_degree_, valency[e] );
    connections_[e] = std::vector< size_t >( valency[e], 0 );
  }
}

//-----------------------------------------------------------------------------

Connectivity::Connectivity(
  std::vector< std::vector< size_t > > const & connectivity )
  : order_( connectivity.size() )
  , min_degree_( ( order_ > 0 ) ? connectivity[0].size() : 0 )
  , max_degree_( 0 )
  , connections_( order_ )
{
  for ( size_t e = 0; e < order_; ++e )
  {
    min_degree_ = std::min( min_degree_, connectivity[e].size() );
    max_degree_ = std::max( max_degree_, connectivity[e].size() );
    connections_[e] = connectivity[e];
  }
}

//-----------------------------------------------------------------------------

Connectivity::Connectivity( Connectivity const & other )
  : order_( other.order_ )
  , min_degree_( other.min_degree_ )
  , max_degree_( other.max_degree_ )
  , connections_( other.order_ )
{
  for ( size_t e = 0; e < order_; ++e )
  {
    connections_[e] = other[e];
  }
}

//-----------------------------------------------------------------------------

Connectivity::~Connectivity() = default;

//-----------------------------------------------------------------------------

auto Connectivity::operator=( Connectivity const & other ) -> Connectivity &
{
  connections_.resize( other.order() );

  for ( size_t e = 0; e < other.order(); ++e )
    connections_[e] = other[e];

  order_      = other.order_;
  min_degree_ = other.min_degree_;
  max_degree_ = other.max_degree_;

#if defined( DEBUG )
  this->check();
#endif

  return *this;
}

//-----------------------------------------------------------------------------

auto Connectivity::operator==( Connectivity const & other ) const -> bool
{
  if ( this == &other )
  {
    return true;
  }
  if ( order_ != other.order_ )
  {
#if DEBUG
    warning( "Connectivity : != order" );
#endif
    return false;
  }
  if ( min_degree_ != other.min_degree_ )
  {
    return false;
  }
  if ( max_degree_ != other.max_degree_ )
  {
    return false;
  }

  if ( connections_.size() != other.connections_.size() )
    return false;

  for ( size_t e = 0; e < order_; ++e )
  {
    if ( connections_[e].size() != other.connections_[e].size() )
      return false;

    for ( size_t f = 0; f < connections_[e].size(); ++f )
    {
      if ( connections_[e][f] != other.connections_[e][f] )
      {
#if DEBUG
        warning( "Connectivity : != connections" );
        message( "First differing entry '[%u][%u]' : %u != %u",
                 e, f, connections_[e][f], other.connections_[e][f] );
#endif
        return false;
      }
    }
  }

  return true;
}


//-----------------------------------------------------------------------------

auto Connectivity::set( std::vector< size_t > const & connectivity ) -> void
{
  if ( connectivity.size() != this->entries() )
  {
    error( "Connectivity : provided connectivity size %u does no match %u",
           connectivity.size(),
           this->entries() );
  }

  size_t pos = 0;
  for ( size_t e = 0; e < order_; ++e )
  {
    std::copy( connectivity.data() + pos,
               connectivity.data() + pos + connections_[e].size(),
               connections_[e].data() );
    pos += connections_[e].size();
  }
}

//-----------------------------------------------------------------------------

auto Connectivity::remap_l( std::vector< size_t > const & map ) -> void
{
  if ( map.size() != order_ )
  {
    error( "Connectivity : remap_left mapping has invalid size" );
  }

  if ( order_ )
  {
    // Set sizes in place of offsets to depict connectivities layout
    std::vector< std::vector< size_t > > remapped( order_ );
    for ( size_t e = 0; e < order_; ++e )
    {
      size_t const ii = map[e];
      dolfin_assert( ii < order_ );
      remapped[e] = connections_[ii];
    }
    std::swap( connections_, remapped );
  }
}

//-----------------------------------------------------------------------------

auto Connectivity::remap_r( std::vector< size_t > const & map ) -> void
{
  for ( size_t e = 0; e < order_; ++e )
  {
    for ( size_t i = 0; i < degree( e ); ++i )
    {
      dolfin_assert( connections_[e][i] < map.size() );
      connections_[e][i] = map[connections_[e][i]];
    }
  }
}

//-----------------------------------------------------------------------------

auto Connectivity::disp() const -> void
{
  section( "Connectivity" );
  message( "order  : %u", order_ );
  message( "degree : %u - %u", min_degree_, max_degree_ );
  end();
}

//-----------------------------------------------------------------------------

auto Connectivity::dump() const -> void
{
  // Display all connections
  for ( size_t e = 0; e < order_; ++e )
  {
    cout << e << ":";
    for ( size_t c = 0; c < connections_[e].size(); ++c )
    {
      cout << " " << connections_[e][c];
    }
    cout << "\n";
  }
}

//-----------------------------------------------------------------------------

auto Connectivity::operator>>( std::vector< size_t > & ) const
  -> Connectivity const &
{
  error( "Connectivity::operator>> unimplemented / deprecated" );
  // A.reserve( this->entries() );
  // for (size_t e = 0; e < order_; ++e)
  //   append( A, connections_[e].begin(), connections_[e].end() );

  // // Set stride if the graph is regular
  // if(min_degree_ == max_degree_)
  //   A %= min_degree_;

  return *this;
}

//-----------------------------------------------------------------------------

auto Connectivity::check() const -> void
{
  /**
   *  CHECK:
   *
   *  Check connectivity size, number of entities and verify that connected
   *  entities are not listed twice.
   *
   */

  for ( size_t e0 = 0; e0 < order_; ++e0 )
  {
    _set< size_t > ce;
    for ( size_t e1 = 0; e1 < connections_[e0].size(); ++e1 )
    {
      size_t ec = connections_[e0][e1];
      if ( ce.count( ec ) > 0 )
      {
        error( "Entity %u appears twice in connectivities for %u", ec, e0 );
      }
      ce.insert( ec );
    }
  }
}

//-----------------------------------------------------------------------------

auto operator<<( std::vector< std::vector< size_t > > & A,
                 Connectivity const &                   C )
  -> std::vector< std::vector< size_t > > &
{
  A.clear();
  A.resize( C.order() );
  for ( size_t e = 0; e < C.order(); ++e )
    A[e] = C[e];
  return A;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
