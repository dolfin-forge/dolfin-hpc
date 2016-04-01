// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-08
// Last changed: 2007-11-30

#ifndef __DOLFIN_MESH_TOPOLOGY_H
#define __DOLFIN_MESH_TOPOLOGY_H

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include "MeshConnectivity.h"
#include "MeshDistributedData.h"
#include "MeshEditor.h"
#include "MeshRenumber.h"
#include "TopologyComputation.h"

namespace dolfin
{

/// MeshTopology stores the topology of a mesh, consisting of mesh entities
/// and connectivity (incidence relations for the mesh entities). Note that
/// the mesh entities do not need to be stored, only the number of entities
/// and the connectivity. Any numbering scheme for the mesh entities is
/// stored separately in a MeshFunction over the entities.
///
/// A mesh entity e may be identified globally as a pair e = (dim, i), where
/// dim is the topological dimension and i is the index of the entity within
/// that topological dimension.

class MeshTopology
{
  friend class TopologyComputation;

public:

  /// Create empty mesh topology
  MeshTopology();

  /// Copy constructor
  MeshTopology(MeshTopology const& topology);

  /// Destructor
  ~MeshTopology();

  /// Assignment
  MeshTopology const& operator=(MeshTopology const& topology);

  /// Equality
  bool operator==(MeshTopology const& other) const;

  /// Non-equality
  bool operator!=(MeshTopology const& other) const;

  /// Return topological dimension
  uint dim() const;

  /// Return number of entities in the local topology for given dimension
  uint size(uint dim) const;

  /// Return number of entities in the global topology for given dimension
  uint global_size(uint dim) const;

  /// Clear all data
  void clear();

  /// Initialize topology of given maximum dimension
  void init(uint dim);

  /// Compute mesh entities of given topological dimension
  uint compute_entities(Mesh& mesh, uint dim) const;

  /// Compute connectivity for given pair of topological dimensions
  void compute_connectivity(Mesh& mesh, uint d0, uint d1) const;

  /// Return connectivity for given pair of topological dimensions
  MeshConnectivity& operator()(uint d0, uint d1);

  /// Return connectivity for given pair of topological dimensions
  MeshConnectivity const& operator()(uint d0, uint d1) const;

  /// Return mesh distribution data
  MeshDistributedData& distdata();

  /// Return mesh distribution data (const version)
  MeshDistributedData const& distdata() const;

  /// Return number of given entities
  uint num_shared(uint dim) const;
  uint num_ghosts(uint dim) const;
  uint num_owned(uint dim) const;

  ///
  void order(Mesh& mesh);

  /// Return if connectivity for given pair is computed
  bool is_computed(uint d0, uint d1) const;

  ///
  bool is_ordered() const;

  ///
  bool is_distributed() const;

  ///
  void renumber(Mesh& mesh);

  /// Return token identifying the internal state of mesh topology
  int token() const;

  /// Display data
  void disp() const;

protected:

  friend void MeshEditor::init_vertices(uint num_vertices);
  friend void MeshEditor::init_cells(uint num_cells);

  /// Set number of entities (size) for given topological dimension
  void init(uint dim, uint size);

private:

  /// Topological dimension
  uint dim_;

  /// Number of mesh entities for each topological dimension
  uint * num_entities_;

  /// Connectivity for pairs of topological dimensions
  MeshConnectivity ** connectivity_;

  /// Distributed mesh topology data
  MeshDistributedData * distdata_;

  /// Return true iff topology is ordered according to the UFC numbering
  bool ordered_;

  //
  int timestamp_;

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
inline MeshConnectivity& MeshTopology::operator()(uint d0, uint d1)
{
  dolfin_assert(d0 <= dim_ && d1 <= dim_);
  return connectivity_[d0][d1];
}

//-----------------------------------------------------------------------------
inline MeshConnectivity const& MeshTopology::operator()(uint d0, uint d1) const
{
  dolfin_assert(d0 <= dim_ && d1 <= dim_);
  return connectivity_[d0][d1];
}

//-----------------------------------------------------------------------------
inline MeshDistributedData& MeshTopology::distdata()
{
  dolfin_assert(distdata_ != NULL);
  return *distdata_;
}

//-----------------------------------------------------------------------------
inline MeshDistributedData const& MeshTopology::distdata() const
{
  dolfin_assert(distdata_ != NULL);
  return *distdata_;
}

//-----------------------------------------------------------------------------
inline uint MeshTopology::global_size(uint dim) const
{
  return (is_distributed() ? distdata().num_global(dim) : this->size(dim));
}

//-----------------------------------------------------------------------------
inline uint MeshTopology::num_shared(uint dim) const
{
  return (is_distributed() ? distdata().num_shared(dim) : 0);
}

//-----------------------------------------------------------------------------
inline uint MeshTopology::num_ghosts(uint dim) const
{
  return (is_distributed() ? distdata().num_ghost(dim) : 0);
}

//-----------------------------------------------------------------------------
inline uint MeshTopology::num_owned(uint dim) const
{
  return (this->size(dim) - this->num_ghosts(dim));
}

//-----------------------------------------------------------------------------
inline bool MeshTopology::is_ordered() const
{
  return ordered_;
}

//-----------------------------------------------------------------------------
inline bool MeshTopology::is_distributed() const
{
  return (distdata_ != NULL);
}

//-----------------------------------------------------------------------------
inline void MeshTopology::renumber(Mesh& mesh)
{
  if (MeshRenumber::renumber(mesh))
  {
    ++renumbering_count_;
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_TOPOLOGY_H */
