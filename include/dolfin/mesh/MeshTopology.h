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
#include "MeshOrdering.h"
#include "MeshRenumber.h"
#include "TopologyComputation.h"

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
  friend class MeshRenumber;
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

  /// Compute mesh entities of given topological dimension
  uint compute_entities(Mesh& mesh, uint dim) const;

  /// Compute connectivity for given pair of topological dimensions
  void compute_connectivity(Mesh& mesh, uint d0, uint d1) const;

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
  void order(Mesh& mesh);

  ///
  bool is_ordered() const;

  ///
  void renumber(Mesh& mesh);

  /// Return token identifying the internal state of mesh topology
  int token() const;

  /// Display data
  void disp() const;

private:

  /// Topological dimension
  uint dim_;

  /// Number of mesh entities for each topological dimension
  uint * num_entities_;

  /// Connectivity for pairs of topological dimensions
  MeshConnectivity ** connectivity_;

  /// Distributed mesh topology data
  MeshDistributedData distdata_;

  /// Return true iff topology is ordered according to the UFC numbering
  bool ordered_;

  //
  int token_;

  //
  uint renumbering_count_;

};

//--- INLINE ------------------------------------------------------------------

//-----------------------------------------------------------------------------
inline uint MeshTopology::dim() const
{
  return dim_;
}

//-----------------------------------------------------------------------------
inline uint MeshTopology::size(uint dim) const
{
  dolfin_assert(dim <= dim_);
  return num_entities_[dim];
}

//-----------------------------------------------------------------------------
inline uint MeshTopology::compute_entities(Mesh& mesh, uint dim) const
{
  return TopologyComputation::computeEntities(mesh, dim);
}

//-----------------------------------------------------------------------------
inline void MeshTopology::compute_connectivity(Mesh& mesh, uint d0, uint d1) const
{
  TopologyComputation::computeConnectivity(mesh, d0, d1);
}

//-----------------------------------------------------------------------------
inline MeshConnectivity& MeshTopology::operator()(uint d0, uint d1)
{
  dolfin_assert(d0 <= dim_ && d1 <= dim_);
  return connectivity_[d0][d1];
}

//-----------------------------------------------------------------------------
const inline MeshConnectivity& MeshTopology::operator()(uint d0, uint d1) const
{
  dolfin_assert(d0 <= dim_ && d1 <= dim_);
  return connectivity_[d0][d1];
}

//-----------------------------------------------------------------------------
inline MeshDistributedData& MeshTopology::distdata()
{
  return distdata_;
}

//-----------------------------------------------------------------------------
const inline MeshDistributedData& MeshTopology::distdata() const
{
  return distdata_;
}

//-----------------------------------------------------------------------------
inline void MeshTopology::order(Mesh& mesh)
{
  MeshOrdering::order(mesh);
  ordered_ = true;
}

//-----------------------------------------------------------------------------
inline bool MeshTopology::is_ordered() const
{
  return ordered_;
}

//-----------------------------------------------------------------------------
inline void MeshTopology::renumber(Mesh& mesh)
{
  if(MeshRenumber::renumber(mesh)) ++renumbering_count_;
}

}

#endif
