// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-08
// Last changed: 2007-11-30

#ifndef __MESH_TOPOLOGY_H
#define __MESH_TOPOLOGY_H

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include "MeshConnectivity.h"
#include "MeshDistributedData.h"

namespace dolfin
{

/// MeshTopology stores the topology of a mesh, consisting of mesh entities
/// and connectivity (incidence relations for the mesh entities). Note that
/// the mesh entities don't need to be stored, only the number of entities
/// and the connectivity. Any numbering scheme for the mesh entities is
/// stored separately in a MeshFunction over the entities.
///
/// A mesh entity e may be identified globally as a pair e = (dim, i), where
/// dim is the topological dimension and i is the index of the entity within
/// that topological dimension.

class MeshTopology
{

  friend class MeshOrdering;
  friend class MPIMeshCommunicator;
  friend class TopologyComputation;

public:

  /// Create empty mesh topology
  MeshTopology();

  /// Copy constructor
  MeshTopology(const MeshTopology& topology);

  /// Destructor
  ~MeshTopology();

  /// Assignment
  const MeshTopology& operator=(const MeshTopology& topology);

  /// Return topological dimension
  uint dim() const;

  /// Return number of entities for given dimension
  uint size(uint dim) const;

  /// Clear all data
  void clear();

  /// Initialize topology of given maximum dimension
  void init(uint dim);

  /// Set number of entities (size) for given topological dimension
  void init(uint dim, uint size);

  /// Return connectivity for given pair of topological dimensions
  MeshConnectivity& operator()(uint d0, uint d1);

  /// Return connectivity for given pair of topological dimensions
  const MeshConnectivity& operator()(uint d0, uint d1) const;

  /// Return mesh distribution data
  MeshDistributedData& distdata();

  /// Return mesh distribution data (const version)
  const MeshDistributedData& distdata() const;

  ///
  bool is_ordered() const;

  /// Display data
  void disp() const;

private:

  /// Topological dimension
  uint _dim;

  /// Number of mesh entities for each topological dimension
  uint * _num_entities;

  /// Connectivity for pairs of topological dimensions
  MeshConnectivity ** _connectivity;

  /// Distributed mesh topology data
  MeshDistributedData _distdata;

  /// Return true iff topology is ordered according to the UFC numbering
  bool _ordered;

};

//--- INLINE ------------------------------------------------------------------

//-----------------------------------------------------------------------------
inline uint MeshTopology::dim() const
{
  return _dim;
}

//-----------------------------------------------------------------------------
inline uint MeshTopology::size(uint dim) const
{
  dolfin_assert(dim <= _dim);
  return _num_entities[dim];
}

//-----------------------------------------------------------------------------
inline MeshConnectivity& MeshTopology::operator()(uint d0, uint d1)
{
  dolfin_assert(d0 <= _dim && d1 <= _dim);
  return _connectivity[d0][d1];
}

//-----------------------------------------------------------------------------
const inline MeshConnectivity& MeshTopology::operator()(uint d0, uint d1) const
{
  dolfin_assert(d0 <= _dim && d1 <= _dim);
  return _connectivity[d0][d1];
}

//-----------------------------------------------------------------------------
inline MeshDistributedData& MeshTopology::distdata()
{
  return _distdata;
}

//-----------------------------------------------------------------------------
const inline MeshDistributedData& MeshTopology::distdata() const
{
  return _distdata;
}

//-----------------------------------------------------------------------------
inline bool MeshTopology::is_ordered() const
{
  return _ordered;
}

}

#endif
