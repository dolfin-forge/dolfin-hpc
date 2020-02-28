// Copyright (C) 2016-2017 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/Connectivity.h>

#include <dolfin/common/Array.h>
#include <dolfin/log/LogStream.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------
Connectivity::Connectivity( uint order, uint degree )
  : order_( order )
  , min_degree_( degree )
  , max_degree_( degree )
  , connections_( order_, Array< uint >( degree, 0 ) )
{
}
//-----------------------------------------------------------------------------
Connectivity::Connectivity( Array< uint > const & valency )
  : order_( valency.size() )
  , min_degree_( ( order_ > 0 ) ? valency[0] : 0 )
  , max_degree_( 0 )
  , connections_( order_ )
{
	for ( uint e = 0; e < order_; ++e )
	{
		min_degree_     = std::min( min_degree_, valency[e] );
		max_degree_     = std::max( max_degree_, valency[e] );
		connections_[e] = Array< uint >( valency[e], 0. );
	}
}
//-----------------------------------------------------------------------------
Connectivity::Connectivity( Array< Array< uint > > const & connectivity )
  : order_( connectivity.size() )
  , min_degree_( ( order_ > 0 ) ? connectivity[0].size() : 0 )
  , max_degree_( 0 )
  , connections_( order_ )
{
	for ( uint e = 0; e < order_; ++e )
	{
		min_degree_ =
		  std::min( min_degree_, static_cast< uint >( connectivity[e].size() ) );
		max_degree_ =
		  std::max( max_degree_, static_cast< uint >( connectivity[e].size() ) );
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
	for ( uint e = 0; e < order_; ++e )
	{
		connections_[e] = other[e];
	}
}
//-----------------------------------------------------------------------------
Connectivity::~Connectivity()
{
}
//-----------------------------------------------------------------------------
Connectivity & Connectivity::operator=( Connectivity const & other )
{
  connections_.resize( other.order() );

  for( uint e = 0; e < other.order(); ++e )
    connections_[e] = other[e];

  order_ = other.order_;
  min_degree_ = other.min_degree_;
  max_degree_ = other.max_degree_;

#if defined( DEBUG )
  this->check();
#endif

  return *this;
}
//-----------------------------------------------------------------------------
bool Connectivity::operator==( Connectivity const & other ) const
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

  for ( uint e = 0; e < order_; ++e )
  {
    if ( connections_[e].size() != other.connections_[e].size() )
      return false;

    for ( uint f = 0; f < connections_[e].size(); ++f )
    {
      if ( connections_[e][f] != other.connections_[e][f] )
      {
#if DEBUG
				warning( "Connectivity : != connections" );
				message( "First differing entry '[%u][%u]' : %u != %u", e, f,
				         connections_[e][f],
				         other.connections_[e][f] );
#endif
        return false;
      }
    }
  }

	return true;
}
//-----------------------------------------------------------------------------
bool Connectivity::operator!=(Connectivity const& other) const
{
  return not (*this == other);
}
//-----------------------------------------------------------------------------
Array< Array< uint > > & Connectivity::operator()()
{
  return connections_;
}
//-----------------------------------------------------------------------------
Array< Array< uint > > const & Connectivity::operator()() const
{
  return connections_;
}
//-----------------------------------------------------------------------------
uint Connectivity::order() const
{
  return order_;
}
//-----------------------------------------------------------------------------
uint Connectivity::entries() const
{
  uint entries = 0;

  for ( uint e = 0; e < order_; ++e )
    entries += connections_[e].size();

  return entries;
}
//----------------------------------------------------------------------------
uint Connectivity::min_degree() const
{
  return min_degree_;
}
//-----------------------------------------------------------------------------
uint Connectivity::max_degree() const
{
  return max_degree_;
}
//-----------------------------------------------------------------------------
uint Connectivity::regular() const
{
  return (min_degree_ == max_degree_ ? min_degree_ : 0);
}
//-----------------------------------------------------------------------------
void Connectivity::set(uint entity, uint const * connections)
{
  dolfin_assert(entity < order_);
  dolfin_assert(not connections_.empty() );
  dolfin_assert( connections != NULL );

  for ( uint e = 0; e < connections_[entity].size(); ++e )
    connections_[entity][e] = connections[e];
}
//-----------------------------------------------------------------------------
void Connectivity::set(Array<uint> const& connectivity)
{
  if (connectivity.size() != this->entries())
  {
    error("Connectivity : provided connectivity size %u does no match %u",
          connectivity.size(), this->entries());
  }

  uint pos = 0;
  for ( uint e = 0; e < order_; ++e )
  {
		std::copy( connectivity.data() + pos,
		           connectivity.data() + pos + connections_[e].size(),
		           connections_[e].data() );
    pos += connections_[e].size();
  }
}
//-----------------------------------------------------------------------------
void Connectivity::remap_l(Array<uint> const& map)
{
  if(map.size() != order_)
  {
    error("Connectivity : remap_left mapping has invalid size");
  }

  if ( order_ )
  {
    // Set sizes in place of offsets to depict connectivities layout
    Array< Array< uint > > remapped( order_ );
    for (uint e = 0; e < order_; ++e)
    {
      uint const ii = map[e];
      dolfin_assert(ii < order_);
      remapped[e] = connections_[ii];
    }
    std::swap(connections_, remapped);
  }
}
//-----------------------------------------------------------------------------
void Connectivity::remap_r(Array<uint> const& map)
{
  for (uint e = 0; e < order_; ++e)
    for ( uint i = 0; i < degree(e); ++i )
      connections_[e][i] = map[connections_[e][i]];
}
//-----------------------------------------------------------------------------
void Connectivity::disp() const
{
  section("Connectivity");
  message("order  : %u", order_);
  message("degree : %u - %u", min_degree_, max_degree_);
  end();
}
//-----------------------------------------------------------------------------
void Connectivity::dump() const
{
  // Display all connections
  for (uint e = 0; e < order_; ++e)
  {
    cout << e << ":";
    for (uint c = 0; c < connections_[e].size(); ++c)
    {
      cout << " " << connections_[e][c];
    }
    cout << "\n";
  }
}
//-----------------------------------------------------------------------------
Connectivity const& Connectivity::operator>>(Array<uint>& A) const
{
  A.reserve( this->entries() );
  for (uint e = 0; e < order_; ++e)
    A.append( connections_[e].begin(), connections_[e].end() );

  // Set stride if the graph is regular
  if(min_degree_ == max_degree_)
    A %= min_degree_;

  return *this;
}
//-----------------------------------------------------------------------------
void Connectivity::check() const
{
  /**
   *  CHECK:
   *
   *  Check connectivity size, number of entities and verify that connected
   *  entities are not listed twice.
   *
   */

  for(uint e0 = 0; e0 < order_; ++e0)
  {
    _set<uint> ce;
    for(uint e1 = 0; e1 < connections_[e0].size(); ++e1)
    {
      uint ec = connections_[e0][e1];
      if(ce.count(ec) > 0)
      {
        error("Entity %u appears twice in connectivities for %u", ec, e0);
      }
      ce.insert(ec);
    }
  }
}
//-----------------------------------------------------------------------------
Array<Array<uint> >& operator<<(Array<Array<uint> >& A, Connectivity const& C)
{
  A.clear();
  A.resize(C.order());
  for (uint e = 0; e < C.order(); ++e)
    A[e] = C[e];
  return A;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

