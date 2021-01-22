// Copyright (C) 2016-2017 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_CONNECTIVITY_H
#define __DOLFIN_MESH_CONNECTIVITY_H

#include <dolfin/common/assert.h>
#include <dolfin/common/types.h>

#include <algorithm>
#include <vector>

namespace dolfin
{

//-----------------------------------------------------------------------------

/**
 *  DOCUMENTATION:
 *
 *  @class  Connectivity
 *
 *  @brief  Connectivity stores a sparse data structure of connections
 *          (incidence relations) between entities.
 *
 */

class Connectivity
{

public:
  /// Create regular connectivity
  Connectivity( size_t order, size_t degree );

  /// Create connectivity
  Connectivity( std::vector< size_t > const & valency );

  /// Create connectivity
  Connectivity( std::vector< std::vector< size_t > > const & connectivity );

  /// Copy constructor
  Connectivity( Connectivity const & other );

  /// Destructor
  ~Connectivity();

  auto operator=( Connectivity const & other ) -> Connectivity &;

  /// Equality
  auto operator==( Connectivity const & other ) const -> bool;

  /// Non-equality
  auto operator!=( Connectivity const & other ) const -> bool;

  /// Return array of connections for given entity
  auto operator[]( size_t entity ) -> std::vector< size_t > &;

  /// Return array of connections for given entity
  auto operator[]( size_t entity ) const -> std::vector< size_t > const &;

  /// Return arrays of connections for all entities
  auto operator()() -> std::vector< std::vector< size_t > > &;

  /// Return arrays of connections for all entities
  auto operator()() const -> std::vector< std::vector< size_t > > const &;

  /// Return incidence of the edge
  auto incident( size_t entity, size_t edge ) const -> bool;

  /// Return index of the edge, -1 if not incident
  auto index( size_t entity, size_t edge ) const -> int;

  //---------------------------------------------------------------------------

  /// Return number of entities
  auto order() const -> size_t;

  /// Return total number of entries
  auto entries() const -> size_t;

  /// Return minimum number of connections
  auto min_degree() const -> size_t;

  /// Return maximum number of connections
  auto max_degree() const -> size_t;

  /// Return degree if regular, zero otherwise
  auto regular() const -> size_t;

  /// Return number of connections for given entity
  auto degree( size_t entity ) const -> size_t;

  /// Set all connections for given entity
  auto set( size_t entity, size_t const * connections ) -> void;

  /// Set all connections for all entities
  auto set( std::vector< size_t > const & connectivity ) -> void;

  /// Remap entities connectivities from old to new ordering, left operator
  auto remap_l( std::vector< size_t > const & map ) -> void;

  /// Remap entities connectivities from old to new ordering, right operator
  auto remap_r( std::vector< size_t > const & map ) -> void;

  /// Display data
  auto disp() const -> void;

  /// Dump data
  auto dump() const -> void;

  //--- SERIALIZATION ---------------------------------------------------------
  auto operator>>( std::vector< size_t > & A ) const -> Connectivity const &;

  //--- CHECK ROUTINES --------------------------------------------------------
  auto check() const -> void;

private:
  /// Number of entities
  size_t order_;

  /// Minimum number of connections
  size_t min_degree_;

  /// Maximum number of connections
  size_t max_degree_;

  /// Connections for all entities stored as a contiguous array
  // size_t **connections_;
  std::vector< std::vector< size_t > > connections_;
};

//-----------------------------------------------------------------------------

inline auto Connectivity::operator!=( Connectivity const & other ) const -> bool
{
  return not( *this == other );
}

//-----------------------------------------------------------------------------

inline auto Connectivity::operator()() -> std::vector< std::vector< size_t > > &
{
  return connections_;
}

//-----------------------------------------------------------------------------

inline auto Connectivity::operator()() const
  -> std::vector< std::vector< size_t > > const &
{
  return connections_;
}

//-----------------------------------------------------------------------------

inline auto Connectivity::order() const -> size_t
{
  return order_;
}

//-----------------------------------------------------------------------------

inline auto Connectivity::entries() const -> size_t
{
  size_t entries = 0;

  for ( size_t e = 0; e < order_; ++e )
    entries += connections_[e].size();

  return entries;
}

//-----------------------------------------------------------------------------

inline auto Connectivity::min_degree() const -> size_t
{
  return min_degree_;
}

//-----------------------------------------------------------------------------

inline auto Connectivity::max_degree() const -> size_t
{
  return max_degree_;
}

//-----------------------------------------------------------------------------

inline auto Connectivity::regular() const -> size_t
{
  return ( min_degree_ == max_degree_ ? min_degree_ : 0 );
}

//-----------------------------------------------------------------------------

inline void Connectivity::set( size_t entity, size_t const * connections )
{
  dolfin_assert( entity < order_ );
  dolfin_assert( not connections_.empty() );
  dolfin_assert( connections != nullptr );

  for ( size_t e = 0; e < connections_[entity].size(); ++e )
    connections_[entity][e] = connections[e];
}

//-----------------------------------------------------------------------------

inline auto Connectivity::degree( size_t entity ) const -> size_t
{
  dolfin_assert( order_ > 0 );
  dolfin_assert( entity < order_ );
  return connections_[entity].size();
}

//-----------------------------------------------------------------------------

inline auto Connectivity::operator[]( size_t entity ) -> std::vector< size_t > &
{
  dolfin_assert( order_ > 0 );
  dolfin_assert( entity < order_ );
  return connections_[entity];
}

//-----------------------------------------------------------------------------

inline auto Connectivity::operator[]( size_t entity ) const
  -> std::vector< size_t > const &
{
  dolfin_assert( order_ > 0 );
  dolfin_assert( entity < order_ );
  return connections_[entity];
}

//-----------------------------------------------------------------------------

inline auto Connectivity::incident( size_t entity, size_t edge ) const -> bool
{
  dolfin_assert( order_ > 0 );
  dolfin_assert( entity < order_ );
  std::vector< size_t >::const_iterator pos =
    std::find( connections_[entity].begin(), connections_[entity].end(), edge );
  return ( connections_[entity].end() != pos );
}

//-----------------------------------------------------------------------------

inline auto Connectivity::index( size_t entity, size_t edge ) const -> int
{
  dolfin_assert( order_ > 0 );
  dolfin_assert( entity < order_ );
  for ( size_t index = 0; index < connections_[entity].size(); ++index )
    if ( connections_[entity][index] == edge )
      return index;
  return -1;
}

//--- OPERATORS ---------------------------------------------------------------

auto operator<<( std::vector< std::vector< size_t > > & A,
                 Connectivity const &                   C )
  -> std::vector< std::vector< size_t > > &;

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_CONNECTIVITY_H */
