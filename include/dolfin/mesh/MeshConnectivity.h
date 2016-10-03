// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-09
// Last changed: 2014-11-03

#ifndef __DOLFIN_MESH_CONNECTIVITY_H
#define __DOLFIN_MESH_CONNECTIVITY_H

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace dolfin
{

template<class T> class Array;

/**
 *  DOCUMENTATION:
 *
 *  @class  MeshConnectivity
 *
 *  @brief  Mesh connectivity stores a sparse data structure of connections
 *          (incidence relations) between mesh entities for a fixed pair of
 *          topological dimensions.
 *          The connectivity can be specified either by first giving the
 *          number of entities and the number of connections for each entity,
 *          which may either be equal for all entities or different, or by
 *          giving the entire (sparse) connectivity pattern.
 *
 */

class MeshConnectivity
{

public:

  /// Create empty connectivity
  MeshConnectivity();

  /// Copy constructor
  MeshConnectivity(MeshConnectivity const& other);

  /// Destructor
  ~MeshConnectivity();

  /// Assignment
  MeshConnectivity const& operator=(MeshConnectivity const& other);

  /// Equality
  bool operator==(MeshConnectivity const& other) const;

  /// Non-equality
  bool operator!=(MeshConnectivity const& other) const;

  /// Initialize number of entities and number of connections (equal for all)
  void init(uint num_entities, uint num_connections);

  /// Initialize number of entities and number of connections (equal for all)
  void init(uint * connectivity, uint num_entities, uint num_connections);

  /// Initialize number of entities and number of connections (individually)
  void init(Array<uint> const& num_connections);

  /// Initialize number of entities and number of connections (individually)
  void init(uint * connectivity, Array<uint> const& num_connections);

  /// Clear all data
  void clear();

  /// Return array of connections for given entity
  uint* operator()(uint entity);

  /// Return array of connections for given entity
  uint const * operator()(uint entity) const;

  /// Return contiguous array of connections for all entities
  uint* operator()();

  /// Return contiguous array of connections for all entities
  uint const * operator()() const;

  //---------------------------------------------------------------------------

  ///
  bool is_initialized() const;

  /// Return number of entities
  uint num_entities() const;

  /// Return total number of connections
  uint size() const;

  /// Return minimum number of connections
  uint min_connections() const;

  /// Return maximum number of connections
  uint max_connections() const;

  /// Return minimum number of connections

  /// Return number of connections for given entity
  uint size(uint entity) const;

  /// Set given connection for given entity
  void set(uint entity, uint connection, uint pos);

  /// Set all connections for given entity
  void set(uint entity, Array<uint> const& connections);

  /// Set all connections for given entity
  void set(uint entity, uint const * connections);

  /// Set all connections for all entities
  void set(Array<uint> const& connectivity);

  /// Set all connections for all entities
  void set(Array<Array<uint> > const& connectivity);

  /// Remap entities connectivities from old to new ordering, left operator
  void remap_left(Array<uint> const& map);

  /// Remap entities connectivities from old to new ordering, right operator
  void remap_right(Array<uint> const& map);

  /// Display data
  void disp() const;

  //--- CHECK ROUTINES --------------------------------------------------------

  /// Check
  void check() const;

private:

  /// Return if initialized
  bool is_initialized_;

  /// Number of entities
  uint num_entities_;

  /// Total number of connections
  uint size_;

  /// Minimum number of connections
  uint min_connections_;

  /// Maximum number of connections
  uint max_connections_;

  /// Offset for first connection for each entity
  uint * offsets_;

  /// Connections for all entities stored as a contiguous array
  uint * connections_;

};

//--- INLINES -----------------------------------------------------------------

inline uint MeshConnectivity::size(uint entity) const
{
  dolfin_assert(num_entities_ > 0);
  dolfin_assert(entity < num_entities_);
  return (offsets_[entity + 1] - offsets_[entity]);
}

//-----------------------------------------------------------------------------
inline uint * MeshConnectivity::operator()(uint entity)
{
  dolfin_assert(num_entities_ > 0);
  dolfin_assert(entity < num_entities_);
  return (connections_ + offsets_[entity]);
}

//-----------------------------------------------------------------------------
inline uint const * MeshConnectivity::operator()(uint entity) const
{
  dolfin_assert(num_entities_ > 0);
  dolfin_assert(entity < num_entities_);
  return (connections_ + offsets_[entity]);
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_CONNECTIVITY_H */
