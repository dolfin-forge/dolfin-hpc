// Copyright (C) 2016-2017 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_CONNECTIVITY_H
#define __DOLFIN_MESH_CONNECTIVITY_H

#include <dolfin/common/types.h>

namespace dolfin
{

template<class T> class Array;

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

  Connectivity & operator=( Connectivity const & other );

  /// Equality
  bool operator==(Connectivity const& other) const;

  /// Non-equality
  bool operator!=(Connectivity const& other) const;

  /// Return array of connections for given entity
  Array< uint > & operator[](uint entity);

  /// Return array of connections for given entity
  Array< uint > const & operator[](uint entity) const;

  /// Return contiguous array of connections for all entities
  Array< Array< uint > > & operator()();

  /// Return contiguous array of connections for all entities
  Array< Array< uint > > const & operator()() const;

  /// Return incidence of the edge
  bool incident(uint entity, uint edge) const;

  /// Return index of the edge, -1 if not incident
  int index(uint entity, uint edge) const;

  //---------------------------------------------------------------------------

  /// Return number of entities
  uint order() const;

  /// Return total number of entries
  uint entries() const;

  /// Return minimum number of connections
  uint min_degree() const;

  /// Return maximum number of connections
  uint max_degree() const;

  /// Return degree if regular, zero otherwise
  uint regular() const;

  /// Return number of connections for given entity
  uint degree(uint entity) const;

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
  Connectivity const& operator>>(Array<uint>& A) const;

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

//--- INLINES -----------------------------------------------------------------

inline uint Connectivity::degree(uint entity) const
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  return connections_[entity].size();
}

//-----------------------------------------------------------------------------
inline Array< uint > & Connectivity::operator[](uint entity)
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  return connections_[entity];
}

//-----------------------------------------------------------------------------
inline Array< uint > const & Connectivity::operator[](uint entity) const
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  return connections_[entity];
}

//-----------------------------------------------------------------------------
inline bool Connectivity::incident(uint entity, uint edge) const
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
	return ( connections_[entity].end()
	         != std::find(
	           connections_[entity].begin(), connections_[entity].end(), edge ) );
}

//-----------------------------------------------------------------------------
inline int Connectivity::index(uint entity, uint edge) const
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  uint index = 0;
  for (; index < connections_[entity].size(); ++index )
    if ( connections_[entity][index] == edge )
      break;

  return ( index == connections_[entity].size() ) ? -1 : index;
}

//--- OPERATORS ---------------------------------------------------------------

Array<Array<uint> >& operator<<(Array<Array<uint> >& A, Connectivity const& C);

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_CONNECTIVITY_H */
