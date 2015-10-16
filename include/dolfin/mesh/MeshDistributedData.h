// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurélien Larcher, 2014.
//
// First added:  2008-07-03
// Last changed: 2014-03-18

#ifndef __DOLFIN_MESH_DISTRIBUTED_DATA_H
#define __DOLFIN_MESH_DISTRIBUTED_DATA_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/EuclideanSpace.h>

namespace dolfin
{

class Vertex;
class Edge;
class Face;
class Mesh;
class MeshEntity;
class MeshTopology;

/**
 *  @class  MeshDistributedData
 *
 *  @brief  Implements data structures for distributed mesh topology.
 */

class MeshDistributedData
{
  /// Set maximum Euclidean dimension
  static uint const MAX_DIMENSION = EuclideanSpace::MAX_DIMENSION;

  /// Set maximum size of data structures
  static uint const MAX_SIZE = MAX_DIMENSION + 1;

public:

  /// Constructor
  MeshDistributedData(MeshTopology& topology);

  /// Destructor
  ~MeshDistributedData();

  /// Assignment
  MeshDistributedData const& operator=(MeshDistributedData const& other);

  /// Return if the distributed data is empty
  bool empty() const;

  /// Return topological dimension of the distributed data
  uint dim() const;

  /// Return if the distributed data is finalized for the given dimension
  bool is_finalized(uint dim) const;

  /// Initialize the distributed data for given topological dimension
  void init(uint dim);

  /// Clear the distributed data
  void clear();

  /// Clear numbering data structures for given dimension (including cache)
  void flush_numbering_data(uint dim);

  /// Clear only numbering C-arrays for given dimension
  void flush_numbering_cache(uint dim);

  /// Clear shared and ghosted data structures for given dimension
  void flush_ownership_data(uint dim);

  /// Clear shared and ghost local-to-adjacent mappings for given dimension
  void flush_mapping(uint dim);

  /// Finalize the data structures for given topological dimension
  void finalize(uint dim);

  //--- Numbering -------------------------------------------------------------

  /// Check if numbering is valid for given topological dimension
  bool has_valid_numbering(uint dim) const;

    /// Invalidate numbering
  void set_invalid_numbering();

  /// Set mapping between local and global numbering
  void set_map(uint local_index, uint global_index, uint dim);

  /// Set new numbering (global-to-local and global-to-local) for given dimension
  void apply_numbering(uint dim, _map<uint, uint> const& local,
                       _map<uint, uint> const& global);

  /// Apply mapping for given dimension
  void apply_mapping(uint dim);

  //--- Global

  /// Return the number of global mesh entities for the given dimension
  uint num_global(uint dim) const;

  /// Return if the global index is registered for given topological dimension
  bool has_global(uint i, uint dim) const;

  /// Return is the given global mesh entity is registered
  bool has_global(MeshEntity const& entity) const;

  /// Return global index of entity given local index and topological dimension
  uint get_global(uint i, uint dim) const;

  /// Return global index of mesh entity
  uint get_global(MeshEntity const& e) const;

  /// Return global index of a vertex given the local index
  uint get_vertex_global(uint i) const;

  /// Return global index of a facet given the local index
  uint get_facet_global(uint i) const;

  /// Return global index of a cell given the local index
  uint get_cell_global(uint i) const;

  /// Set the number of global entities for given topological dimension
  /// TODO: Deprecate.
  void set_num_global(uint dim, uint const num_global);

  /// Compute the number of global entities for given topological dimension
  void apply_num_global(uint dim, uint& offset);

  //--- Local

  /// Return the number of local mesh entities for the given dimension
  uint num_local(uint dim) const;

  /// Return if the local index is registered for given topological dimension
  bool has_local(uint i, uint dim) const;

  /// Return if the given local mesh entity is registered
  bool has_local(MeshEntity const& entity) const;

  /// Return local index of entity given global index and topological dimension
  uint get_local(uint i, uint dim) const;

  /// Return local index of mesh entity
  uint get_local(MeshEntity const& e) const;

  /// Return local index of a vertex given the global index
  uint get_vertex_local(uint i) const;

  /// Return local index of a facet given the global index
  uint get_facet_local(uint i) const;

  /// Return local index of a cell given the global index
  uint get_cell_local(uint i) const;

  //--- Ownership -------------------------------------------------------------

  /// Check if numbering is valid for given topological dimension
  bool has_valid_ownership(uint dim) const;

  /// Check if mapping is valid for given topological dimension
  bool has_valid_mapping(uint dim) const;

  /// Invalidate ownership data structures
  void set_invalid_ownership();

  /// Apply ownership data (currently just set to valid)
  void apply_ownership(uint dim);

  /// Remap ghosted mesh entities' owner identities
  void remap_ownership(int const* mapping);

  //--- Sharedness ---

  /// Return if the mesh entity for given index and dimension is shared
  bool is_shared(uint i, uint dim) const;

  /// Return if the given mesh entity is shared
  bool is_shared(MeshEntity const& entity) const;

  /// Return the number of shared entities for the given topological dimension
  uint num_shared(uint dim) const;

  /// Return the number of shared facets per adjacent rank
  uint num_shared_with(uint rank, uint dim) const;

  /// Return the mapping of shared entities ordering from local ordering of
  /// shared iterator to adjacent rank ordering (send)
  Array<uint> const& get_shared_mapping_to(uint rank, uint dim) const;
  Array<uint>& get_shared_mapping_to(uint rank, uint dim);

  /// Return the mapping of shared entities ordering from adjacent rank ordering
  /// local ordering of shared iterator (recv)
  Array<uint> const& get_shared_mapping_from(uint rank, uint dim) const;
  Array<uint>& get_shared_mapping_from(uint rank, uint dim);

  /// Set mesh entity for given local index and topological dimension as shared
  void set_shared(uint local_index, uint dim);

  /// Set the given mesh entity as shared
  void set_shared(MeshEntity const& m);

  //--- Adjacency ---

  /// Return the set of adjacent ranks for the given topological dimension
  _set<uint> const& get_adj_ranks(uint dim) const;

  /// Return the number of adjacent ranks for the given dimension
  uint num_adj_ranks(uint dim) const;

  /// Return the set of adjacent ranks of a shared entity given its local index
  /// and topological dimension
  _set<uint> const& get_shared_adj(uint local_index, uint dim) const;

  /// Return the set of adjacent ranks of a shared entity
  _set<uint> const& get_shared_adj(MeshEntity const& m) const;

  /// Set rank as adjacent for the mesh entity given index and dimension
  void set_shared_adj(uint i, uint rank, uint dim);

  /// Set rank as adjacent for the given mesh entity
  void set_shared_adj(MeshEntity const& m, uint rank);

  /// Set all adjacent ranks for the mesh entity given index and dimension
  /// Previous adjacency is erased
  void setall_shared_adj(uint i, _set<uint> const& ranks, uint dim);

  /// Set all adjacent ranks for the given mesh entity
  /// Previous adjacency is erased
  void setall_shared_adj(MeshEntity const& m, _set<uint> const& ranks);

  //--- Ghostedness ---

  /// Return if the mesh entity for given index and dimension is ghosted
  bool is_ghost(uint i, uint dim) const;

  /// Return if the given mesh entity is ghosted
  bool is_ghost(MeshEntity const& entity) const;

  /// Return the number of owned entities for the given topological dimension
  uint num_owned(uint dim) const;

  /// Return the number of ghosted entities for the given topological dimension
  uint num_ghost(uint dim) const;

  /// Return the number of ghosted entities per owner rank
  uint num_ghost_from(uint rank, uint dim) const;

  /// Return owner of the ghosted entity given its global index and topological
  /// dimension
  uint get_owner(uint local_index, uint dim) const;

  /// Return owner of the ghosted entity
  uint get_owner(MeshEntity const& m) const;

  /// Return the mapping of ghost entities ordering from local ordering of
  /// ghost iterator to adjacent rank ordering (send)
  Array<uint> const& get_ghost_mapping_to(uint rank, uint dim) const;
  Array<uint>& get_ghost_mapping_to(uint rank, uint dim);

  /// Return the mapping of ghost entities ordering from adjacent rank ordering
  /// local ordering of ghost iterator (recv)
  Array<uint> const& get_ghost_mapping_from(uint rank, uint dim) const;
  Array<uint>& get_ghost_mapping_from(uint rank, uint dim);

  /// Set mesh entity for given local index and topological dimension as ghosted
  void set_ghost(uint local_index, uint dim);

  /// Set the given mesh entity as ghosted
  void set_ghost(MeshEntity const& m);

  /// Set the owner of the ghosted mesh entity with given local index and
  /// topological dimension
  void set_ghost_owner(uint i, uint rank, uint dim);

  /// Set the owner of the given ghosted mesh entity
  void set_ghost_owner(MeshEntity const& m, uint rank);

  //---

  /// Display basic information
  void disp() const;

  //--- Debugging

  ///
  bool check(bool throw_error = false) const;

  ///
  bool check_shared(uint local_index, uint dim, bool error = false) const;

  ///
  bool check_ghost(uint local_index, uint dim, bool error = false) const;

protected:

private:

  friend class MeshGhostIterator;
  friend class MeshSharedIterator;

  MeshTopology& topology_;

  // Topological dimensions
  uint topological_dim_;
  mutable uint cell_dim_;
  mutable uint facet_dim_;

  // Numbering for entities of topological dimension
  bool valid_numbering_[MAX_SIZE];
  uint num_global_[MAX_SIZE];
  mutable _map<uint, uint> global_indices_[MAX_SIZE];
  mutable _map<uint, uint> local_indices_[MAX_SIZE];
  uint max_global_index_;

  // Array caching of numbering
  bool finalized_[MAX_SIZE];
  uint cached_global_size_[MAX_SIZE];
  uint * cached_global_indices_[MAX_SIZE];

  // Ownership
  // NOTE: Cells are not shared/ghosted for the moment but the data structure
  //       is in place.
  bool valid_ownership_[MAX_SIZE];
  _set<uint> shared_[MAX_SIZE];
  _set<uint> adjacent_ranks_[MAX_SIZE];
  mutable _map<uint, _set<uint> > shared_adj_[MAX_SIZE];
  _set<uint> ghost_[MAX_SIZE];
  mutable _map<uint, uint> ghost_owner_[MAX_SIZE];

  // Adjacent mappings and reverse mappings
  bool valid_mapping_[MAX_SIZE];
  typedef _map<uint, std::pair< Array<uint>, Array<uint> > > AdjacentMapping;
  AdjacentMapping shared_mapping_[MAX_SIZE];
  AdjacentMapping ghost_mapping_[MAX_SIZE];

};

//-----------------------------------------------------------------------------

/**
 *  @class  MeshGhostIterator
 *
 *  @brief  Implements an iterator on ghosted mesh entities.
 */

class MeshGhostIterator
{

public:

  MeshGhostIterator(MeshDistributedData const& distdata, uint i) :
      distdata_(distdata),
      dim_(i)
  {
    iter_ = distdata_.ghost_[i].begin();
  }

  ~MeshGhostIterator()
  {
  }

  MeshGhostIterator& operator++()
  {
    ++iter_;
    return *this;
  }

  inline uint index() const
  {
    return *iter_;
  }

  inline uint owner() const
  {
    return distdata_.get_owner(*iter_, dim_);
  }

  inline bool end() const
  {
    return iter_ == distdata_.ghost_[dim_].end();
  }

  inline _set<uint> const& adj() const
  {
    return distdata_.shared_adj_[dim_][*iter_];
  }

  //--- Debugging ---

  inline bool check() const
  {
    return distdata_.check_ghost(*iter_, dim_);
  }

private:

  MeshDistributedData const& distdata_;
  uint const dim_;
  _set<uint>::const_iterator iter_;

};

//-----------------------------------------------------------------------------

/**
 *  @class  MeshSharedIterator
 *
 *  @brief  Implements an iterator on shared mesh entities.
 */

class MeshSharedIterator
{

public:

  MeshSharedIterator(MeshDistributedData const& distdata, uint i) :
      distdata_(distdata),
      dim_(i)
  {
    dolfin_assert(i <= distdata.dim());
    iter_ = distdata_.shared_[i].begin();
  }

  ~MeshSharedIterator()
  {
  }

  MeshSharedIterator& operator++()
  {
    ++iter_;
    return *this;
  }

  inline uint index() const
  {
    return *iter_;
  }

  inline bool end() const
  {
    return iter_ == distdata_.shared_[dim_].end();
  }

  inline _set<uint> const& adj() const
  {
    return distdata_.shared_adj_[dim_][*iter_];
  }

  //--- Debugging ---

  inline bool check() const
  {
    return distdata_.check_shared(*iter_, dim_);
  }

private:

  MeshDistributedData const& distdata_;
  uint const dim_;
  _set<uint>::const_iterator iter_;
};

}

#endif
