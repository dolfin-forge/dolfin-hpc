// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurélien Larcher, 2014.
//
// First added:  2008-07-03
// Last changed: 2014-03-18

#ifndef __MESH_DISTRIBUTED_DATA_H
#define __MESH_DISTRIBUTED_DATA_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

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
 *  @brief
 */

class MeshDistributedData
{
  static uint const MAX_DIM = 3;

public:

  /// Constructor
  MeshDistributedData(MeshTopology& topology);

  /// Destructor
  ~MeshDistributedData();

  /// Assignment
  MeshDistributedData const& operator=(MeshDistributedData const& other);

  /// Return if the distributed data is empty
  bool empty() const;

  /// Initialize the distributed data for given topological dimension
  void init(uint const dim);

  /// Clear the distributed data
  void clear();

  /// Finalize the data structures for given topological dimension
  void finalize(uint const dim);

  //--- Numbering -------------------------------------------------------------

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

  /// Return the number of global entities for the given dimension
  uint num_global(uint dim) const;

  /// Set mapping between local and global nun
  void set_map(uint local_index, uint global_index, uint dim);

  /// Set the number of global entities for given topological dimension
  void set_num_global(uint dim, uint num_global);

  /// Invalidate numbering
  void set_invalid_numbering();

  //--- Ownership -------------------------------------------------------------

  //--- Adjacency ---

  /// Return the set of adjacent ranks for the given topological dimension
  _set<uint> const& get_adj(uint dim) const;

  /// Return the number of adjacent ranks for the given dimension
  uint num_adj(uint dim) const;

  //--- Sharedness ---

  ///
  bool is_shared(uint i, uint dim) const;

  ///
  bool is_shared(MeshEntity const& entity) const;

  ///
  uint num_shared(uint dim) const;

  /// Return the number of shared facets per adjacent rank
  uint num_shared_with(uint rank, uint dim) const;

  /// Return the set of adjacent ranks of a shared entity given its local index
  /// and topological dimension
  _set<uint> const& get_shared_adj(uint local_index, uint dim) const;

  /// Return the set of adjacent ranks of a shared entity
  _set<uint> const& get_shared_adj(MeshEntity const& m) const;

  /// Return the mapping of shared entities ordering from local ordering of
  /// shared iterator to adjacent rank ordering (send), and the converse (recv)
  Array<uint> const& get_shared_mapping_to(uint rank, uint dim) const;
  Array<uint> const& get_shared_mapping_from(uint rank, uint dim) const;

  void set_shared(uint local_index, uint dim);
  void set_shared(MeshEntity const& m);

  void set_shared_adj(uint i, uint rank, uint dim);
  void set_shared_adj(MeshEntity const& m, uint rank);

  void setall_shared_adj(uint i, _set<uint> const& ranks, uint dim);
  void setall_shared_adj(MeshEntity const& m, _set<uint> const& ranks);

  //--- Ghostedness ---

  inline bool is_ghost(uint i, uint dim) const
  { return (MPI::numProcesses() > 1 ? (ghost_[dim].count(i) > 0) : false);}

  bool is_ghost(MeshEntity const& entity) const;

  inline uint num_ghost(uint dim) const
  { return ghost_[dim].size();}

  /// Return owner of the ghosted entity given its global index and topological
  /// dimension
  uint get_owner(uint local_index, uint dim) const;

  /// Return owner of the ghosted entity
  uint get_owner(MeshEntity const& m) const;

  /// Return the mapping of ghost_ entities ordering from local ordering of
  /// shared iterator to adjacent rank ordering (send), and the converse (recv)
  Array<uint> const& get_ghost_mapping_to(uint rank, uint dim) const;
  Array<uint> const& get_ghost_mapping_from(uint rank, uint dim) const;

  /// Return the number of ghosted facets per owner rank
  uint num_ghost_from(uint rank, uint dim) const;

  //---

  void set_ghost(uint local_index, uint dim);
  void set_ghost(MeshEntity const& m);

  void set_ghost_owner(uint i, uint rank, uint dim);
  void set_ghost_owner(MeshEntity const& m, uint rank);

  inline void set_invalid_ownership()
  {
    set_invalid_numbering();
    flush_edges(); flush_faces();
    finalized_ = false;
  }

  void remap_owner(int* mapping);

  inline void flush_edges()
  { shared_[1].clear(); ghost_[1].clear(); ghost_owner_[1].clear();}

  inline void flush_faces()
  { shared_[2].clear(); ghost_[2].clear(); ghost_owner_[2].clear();}

  void flush_mappings(uint dim);

  //
  void disp() const;

protected:

private:

  uint topological_dim_;
  mutable uint cell_dim_;
  mutable uint facet_dim_;

  uint max_global_vertex_index_;
  uint num_global_[MAX_DIM+1];
  bool valid_numbering_[MAX_DIM+1];

  bool valid_edge_ownership_;
  bool valid_face_ownership_;

  bool valid_shared_facets_mapping_;

  mutable _map<uint, uint> global_indices_[MAX_DIM+1];
  mutable _map<uint, uint> local_indices_[MAX_DIM+1];

  mutable _map<uint, uint> ghost_owner_[MAX_DIM];
  mutable _map<uint, _set<uint> > shared_adj_[MAX_DIM];

  _set<uint> adjacent_ranks_[MAX_DIM];
  _set<uint> shared_[MAX_DIM];
  _set<uint> ghost_[MAX_DIM];

  // Adjacent mappings and reverse mappings
  typedef _map<uint, std::pair< Array<uint>, Array<uint> > > AdjacentMapping;
  AdjacentMapping shared_mapping_[MAX_DIM];
  AdjacentMapping ghost_mapping_[MAX_DIM];

  bool finalized_;
  uint * global_vertex_indices_;
  uint * global_facet_indices_;
  uint * global_cell_indices_;

  uint global_vertex_indices_size_;
  uint global_facet_indices_size_;
  uint global_cell_indices_size_;

  friend class MeshGhostIterator;
  friend class MeshSharedIterator;
  friend class MeshRenumber;

};

class MeshGhostIterator
{

public:

  MeshGhostIterator(MeshDistributedData& distdata, uint i) :
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

private:

  MeshDistributedData& distdata_;
  uint const dim_;
  _set<uint>::iterator iter_;

};

class MeshSharedIterator
{

public:

  MeshSharedIterator(MeshDistributedData& distdata, uint i) :
      distdata_(distdata),
      dim_(i)
  {
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

  inline _set<uint> adj() const
  {
    return distdata_.shared_adj_[dim_][*iter_];
  }

private:

  MeshDistributedData& distdata_;
  uint const dim_;
  _set<uint>::iterator iter_;
};

}

#endif
