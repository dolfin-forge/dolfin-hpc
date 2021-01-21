// Copyright (C) 2016-2017 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_CONNECTIVITY_H
#define __DOLFIN_MESH_CONNECTIVITY_H

#include <dolfin/common/assert.h>
#include <dolfin/common/Array.h>
#include <dolfin/common/types.h>


namespace dolfin
{

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
  Connectivity(uint order, uint degree);

  /// Create connectivity
  Connectivity(Array<uint> const& valency);

  /// Create connectivity
  Connectivity(Array< Array<uint> > const& connectivity);

  /// Copy constructor
  Connectivity(Connectivity const& other);

  /// Destructor
  ~Connectivity();

  auto operator=( Connectivity const & other ) -> Connectivity &;

  /// Equality
  auto operator==(Connectivity const& other) const -> bool;

  /// Non-equality
  auto operator!=(Connectivity const& other) const -> bool;

  /// Return array of connections for given entity
  auto operator[](uint entity) -> Array< uint > &;

  /// Return array of connections for given entity
  auto operator[](uint entity) const -> Array< uint > const &;

  /// Return arrays of connections for all entities
  auto operator()() -> Array< Array< uint > > &;

  /// Return arrays of connections for all entities
  auto operator()() const -> Array< Array< uint > > const &;

  /// Return incidence of the edge
  auto incident(uint entity, uint edge) const -> bool;

  /// Return index of the edge, -1 if not incident
  auto index(uint entity, uint edge) const -> int;

  //---------------------------------------------------------------------------

  /// Return number of entities
  auto order() const -> uint;

  /// Return total number of entries
  auto entries() const -> uint;

  /// Return minimum number of connections
  auto min_degree() const -> uint;

  /// Return maximum number of connections
  auto max_degree() const -> uint;

  /// Return degree if regular, zero otherwise
  auto regular() const -> uint;

  /// Return number of connections for given entity
  auto degree(uint entity) const -> uint;

  /// Set all connections for given entity
  void set(uint entity, uint const * connections);

  /// Set all connections for all entities
  void set(Array<uint> const& connectivity);

  /// Remap entities connectivities from old to new ordering, left operator
  void remap_l(Array<uint> const& map);

  /// Remap entities connectivities from old to new ordering, right operator
  void remap_r(Array<uint> const& map);

  /// Display data
  void disp() const;

  /// Dump data
  void dump() const;

  //--- SERIALIZATION ---------------------------------------------------------
  auto operator>>(Array<uint>& A) const -> Connectivity const&;

  //--- CHECK ROUTINES --------------------------------------------------------
  void check() const;

private:
  /// Number of entities
  uint order_;

  /// Minimum number of connections
  uint min_degree_;

  /// Maximum number of connections
  uint max_degree_;

  /// Connections for all entities stored as a contiguous array
  // uint **connections_;
  Array< Array< uint > > connections_;

};

//-----------------------------------------------------------------------------
inline auto Connectivity::operator!=(Connectivity const& other) const -> bool
{
  return not (*this == other);
}

//-----------------------------------------------------------------------------
inline auto Connectivity::operator()() -> Array< Array< uint > > &
{
  return connections_;
}

//-----------------------------------------------------------------------------
inline auto Connectivity::operator()() const -> Array< Array< uint > > const &
{
  return connections_;
}

//-----------------------------------------------------------------------------
inline auto Connectivity::order() const -> uint
{
  return order_;
}

//-----------------------------------------------------------------------------
inline auto Connectivity::entries() const -> uint
{
  uint entries = 0;

  for ( uint e = 0; e < order_; ++e )
    entries += connections_[e].size();

  return entries;
}

//----------------------------------------------------------------------------
inline auto Connectivity::min_degree() const -> uint
{
  return min_degree_;
}

//-----------------------------------------------------------------------------
inline auto Connectivity::max_degree() const -> uint
{
  return max_degree_;
}

//-----------------------------------------------------------------------------
inline auto Connectivity::regular() const -> uint
{
  return (min_degree_ == max_degree_ ? min_degree_ : 0);
}

//-----------------------------------------------------------------------------
inline void Connectivity::set(uint entity, uint const * connections)
{
  dolfin_assert(entity < order_);
  dolfin_assert(not connections_.empty() );
  dolfin_assert( connections != nullptr );

  for ( uint e = 0; e < connections_[entity].size(); ++e )
    connections_[entity][e] = connections[e];
}

inline auto Connectivity::degree(uint entity) const -> uint
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  return connections_[entity].size();
}

//-----------------------------------------------------------------------------
inline auto Connectivity::operator[](uint entity) -> Array< uint > &
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  return connections_[entity];
}

//-----------------------------------------------------------------------------
inline auto Connectivity::operator[](uint entity) const -> Array< uint > const &
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  return connections_[entity];
}

//-----------------------------------------------------------------------------
inline auto Connectivity::incident(uint entity, uint edge) const -> bool
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  Array< uint >::const_iterator pos =
    std::find( connections_[entity].begin(), connections_[entity].end(), edge );
  return ( connections_[entity].end() != pos );
}

//-----------------------------------------------------------------------------
inline auto Connectivity::index(uint entity, uint edge) const -> int
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  for (uint index = 0; index < connections_[entity].size(); ++index )
    if ( connections_[entity][index] == edge )
      return index;
  return -1;
}

//--- OPERATORS ---------------------------------------------------------------

auto operator<<(Array<Array<uint> >& A, Connectivity const& C) -> Array<Array<uint> >&;

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_CONNECTIVITY_H */
