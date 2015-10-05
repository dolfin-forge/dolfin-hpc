// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-09
// Last changed: 2014-11-03

#ifndef __DOLFIN_MESH_CONNECTIVITY_H
#define __DOLFIN_MESH_CONNECTIVITY_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>

namespace dolfin
{

/// Mesh connectivity stores a sparse data structure of connections
/// (incidence relations) between mesh entities for a fixed pair of
/// topological dimensions.
///
/// The connectivity can be specified either by first giving the
/// number of entities and the number of connections for each entity,
/// which may either be equal for all entities or different, or by
/// giving the entire (sparse) connectivity pattern.

class MeshConnectivity
{
public:

  /// Create empty connectivity
  MeshConnectivity();

  /// Copy constructor
  MeshConnectivity(const MeshConnectivity& other);

  /// Destructor
  ~MeshConnectivity();

  /// Assignment
  MeshConnectivity const& operator=(MeshConnectivity const& other);

  /// Return total number of connections
  uint size() const;

  /// Return number of connections for given entity
  uint size(uint entity) const;

  /// Return array of connections for given entity
  uint* operator()(uint entity);

  /// Return array of connections for given entity
  uint const * operator()(uint entity) const;

  /// Return contiguous array of connections for all entities
  uint* operator()();

  /// Return contiguous array of connections for all entities
  uint const * operator()() const;

  /// Clear all data
  void clear();

  /// Initialize number of entities and number of connections (equal for all)
  void init(uint num_entities, uint num_connections);

  /// Initialize number of entities and number of connections (individually)
  void init(Array<uint> const& num_connections);

  /// Set given connection for given entity
  void set(uint entity, uint connection, uint pos);

  /// Set all connections for given entity
  void set(uint entity, Array<uint> const& connections);

  /// Set all connections for given entity
  void set(uint entity, uint const * connections);

  /// Set all connections for all entities
  void set(Array<Array<uint> > const& connectivity);

  /// Display data
  void disp() const;

private:

  friend class MPIMeshCommunicator;

  /// Total number of connections
  uint size_;

  /// Number of entities
  uint num_entities_;

  /// Connections for all entities stored as a contiguous array
  uint * connections_;

  /// Offset for first connection for each entity
  uint * offsets_;

};

//--- INLINES -----------------------------------------------------------------

inline uint MeshConnectivity::size() const
{
  return size_;
}

//-----------------------------------------------------------------------------
inline uint MeshConnectivity::size(uint entity) const
{
  return (entity < num_entities_ ? offsets_[entity + 1] - offsets_[entity] : 0);
}

//-----------------------------------------------------------------------------
inline uint * MeshConnectivity::operator()(uint entity)
{
  return (entity < num_entities_ ? connections_ + offsets_[entity] : 0);
}

//-----------------------------------------------------------------------------
inline const uint * MeshConnectivity::operator()(uint entity) const
{
  return (entity < num_entities_ ? connections_ + offsets_[entity] : 0);
}

//-----------------------------------------------------------------------------
inline uint * MeshConnectivity::operator()()
{
  return connections_;
}

//-----------------------------------------------------------------------------
inline const uint * MeshConnectivity::operator()() const
{
  return connections_;
}

//-----------------------------------------------------------------------------

}

#endif
