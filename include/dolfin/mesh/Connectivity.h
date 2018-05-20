// Copyright (C) 2016-2017 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __DOLFIN_MESH_CONNECTIVITY_H
#define __DOLFIN_MESH_CONNECTIVITY_H

#include <dolfin/common/Clonable.h>

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

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

class Connectivity : public Clonable<Connectivity>
{

public:

  /// Create regular connectivity
  Connectivity(uint order, uint degree);

  /// Create connectivity
  Connectivity(Array<uint> const& valency);

  /// Create connectivity
  Connectivity(Array<Array<uint> > const& connectivity);

  /// Copy constructor
  Connectivity(Connectivity const& other);

  /// Destructor
  ~Connectivity();

  /// Equality
  bool operator==(Connectivity const& other) const;

  /// Non-equality
  bool operator!=(Connectivity const& other) const;

  /// Return array of connections for given entity
  uint* operator()(uint entity);

  /// Return array of connections for given entity
  uint const * operator()(uint entity) const;

  /// Set array bounds of connections for given entity
  void operator()(uint entity, uint *& b, uint *& e);

  /// Set array bounds of connections for given entity
  void operator()(uint entity, uint const *& b, uint const *& e) const;

  /// Return contiguous array of connections for all entities
  uint* operator()();

  /// Return contiguous array of connections for all entities
  uint const * operator()() const;

  /// Return incidence of the edge
  bool incident(uint entity, uint edge) const;

  /// Return index of the edge, -1 if not incident
  int index(uint entity, uint edge) const;

  //---------------------------------------------------------------------------

  /// Return number of entities
  uint order() const;

  /// Return total number of entries
  uidx entries() const;

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

  //--- ITERATOR --------------------------------------------------------------

  typedef uint*       data_iterator;
  typedef uint const* const_data_iterator;

  inline data_iterator data()  const { return connections_[0]; }
  inline data_iterator bound() const { return connections_[order_]; }

  //--- SERIALIZATION ---------------------------------------------------------
  Connectivity const& operator>>(Array<uint>& A) const;

  //--- CHECK ROUTINES --------------------------------------------------------

  /// Check
  void check() const;

private:

  /// Create empty connectivity (Disabled)
  Connectivity();

  /// Number of entities
  uint order_;

  /// Minimum number of connections
  uint min_degree_;

  /// Maximum number of connections
  uint max_degree_;

  /// Connections for all entities stored as a contiguous array
  uint **connections_;

};

//--- INLINES -----------------------------------------------------------------

inline uint Connectivity::degree(uint entity) const
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  return (connections_[entity + 1] - connections_[entity]);
}

//-----------------------------------------------------------------------------
inline uint * Connectivity::operator()(uint entity)
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  return connections_[entity];
}

//-----------------------------------------------------------------------------
inline uint const * Connectivity::operator()(uint entity) const
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
  uint const * e = connections_[entity];
  uint const * const n = connections_[entity + 1];
  while (e != n && *e != edge) { ++e; }
  return (e != n);
}

//-----------------------------------------------------------------------------
inline int Connectivity::index(uint entity, uint edge) const
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  uint const * const b = connections_[entity];
  uint const * const n = connections_[entity + 1];
  uint const * e = b;
  while (e != n && *e != edge) { ++e; }
  return (e == n ? -1 : (e - b));
}

//-----------------------------------------------------------------------------
inline void Connectivity::operator()(uint entity, uint *& b, uint *& e)
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  b = connections_[entity];
  e = connections_[entity + 1];
}

//-----------------------------------------------------------------------------
inline void Connectivity::operator()(uint entity, uint const *& b,
                                         uint const *& e) const
{
  dolfin_assert(order_ > 0);
  dolfin_assert(entity < order_);
  b = connections_[entity];
  e = connections_[entity + 1];
}

//--- OPERATORS ---------------------------------------------------------------

Array<Array<uint> >& operator<<(Array<Array<uint> >& A, Connectivity const& C);

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_CONNECTIVITY_H */
