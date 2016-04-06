//
//
//

#ifndef __DOLFIN_DISTRIBUTED_DATA_H
#define __DOLFIN_DISTRIBUTED_DATA_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

/*
 *  @class  DistributedData
 */

class DistributedData
{

  friend class SharedIterator;
  friend class GhostIterator;

public:

  ///
  bool valid_numbering;
  bool valid_ownership;
  bool valid_adjacency;

  ///
  DistributedData();

  ///
  DistributedData(DistributedData const& other);

  ///
  ~DistributedData();

  ///
  DistributedData& operator=(DistributedData const& other);

  ///
  bool operator==(DistributedData const& other) const;

  ///
  bool operator!=(DistributedData const& other) const;

  /// Clear all data
  void clear();

  /// Finalize the data: validate and set process range + global size
  void finalize();

  /// Return the storage capacity of mappings
  uint capacity() const;

  /// Return whether the data is finalized
  bool is_finalized() const;

  //--- Data bounds -----------------------------------------------------------

  /// Return the offset of the process range
  uint offset() const;

  /// Return the local data size
  uint range_size() const;

  /// Return if the global index is in the process range
  bool in_range(uint global_index) const;

  /// Return the local data size
  uint local_size() const;

  /// Return the global data size i.e the global sum of number of owned entities
  uint global_size() const;

  //--- Numbering -------------------------------------------------------------

  /// Return if the index is a local index
  uint has_local(uint local_index) const;

  /// Return the global index given a local index
  uint get_global(uint local_index) const;

  /// Return if the index is a global index
  uint has_global(uint global_index) const;

  /// Return the local index given a global index
  uint get_local(uint global_index) const;

  /// Set local-to-global mapping
  void set_map(uint local_index, uint global_index);

  /// Re-map numbering with given mapping
  void remap_numbering(Array<uint> const& mapping);

  /// Re-index global indices to have contiguous numbering of owned entities
  void renumber_global();

  //--- Adjacency -------------------------------------------------------------

  /// Return whether the given rank is adjacent
  bool has_adj_rank(uint rank) const;

  /// Return the number of adjacent ranks
  uint num_adj_ranks() const;

  /// Return the set of adjacent ranks
  _set<uint> const& get_adj_ranks() const;

  //--- Ownership -------------------------------------------------------------

  /// Return the owner of the entity: self if owned and not self otherwise
  uint get_owner(uint local_index) const;

  /// Return if the entity is owned
  bool is_owned(uint local_index) const;

  /// Return if the entity is shared: it can be owned or not
  bool is_shared(uint local_index) const;

  /// Return if the entity is ghost:  it is shared and not owned
  bool is_ghost(uint local_index) const;

  /// Return the number of owned entities
  uint num_owned() const;

  /// Return the number of shared entities
  uint num_shared() const;

  /// Return the number of ghost entities
  uint num_ghost() const;

  /// Re-map ownership with given mapping
  void remap_ownership(Array<uint> const& mapping);

  //--- Shared ---

  /// Return the adjacent set of a shared entity
  _set<uint> const& get_shared_adj(uint local_index) const;

  /// Return the common adjacent set to an array of shared entities
  void get_common_adj(uint n, uint const * indices, _set<uint>& adjs) const;

  /// Set the entity as shared, the adjacent set is not created.
  /// If the entity is ghost then it stays that way
  void set_shared(uint local_index);

  /// Add a rank as adjacent, this cannot be self
  void set_shared_adj(uint local_index, uint adj);

  /// Set the adjacent set for the given shared entity, this cannot contain self
  void setall_shared_adj(uint local_index, _set<uint> const& adjs);

  //--- Ghosts ---

  /// Set the given entity as ghost
  void set_ghost(uint local_index, uint owner);

  //
  void disp() const;

private:

  uint rank_;
  uint pe_size_;
  uint offset_;
  uint range_size_;
  uint global_size_;

  //
  bool finalized_;

  ///
  _map<uint, uint> global_;
  _map<uint, uint> local_;

  //
  _set<uint> adjacents_;
  _map<uint, _set<uint> > shared_;
  _map<uint, uint> ghost_;

  ///
  uint * cached_numbering_;
  uint * cached_ownership_;

};

/**
 *  @class  SharedIterator
 *
 *  @brief  Implements an iterator on shared entities.
 */

class SharedIterator
{

public:

  ///
  SharedIterator(DistributedData const& distdata) :
      distdata_(distdata),
      iter_(distdata_.shared_.begin())
  {
  }

  ///
  ~SharedIterator()
  {
  }

  ///
  SharedIterator& operator++()
  {
    ++iter_;
    return *this;
  }

  ///
  inline uint index() const
  {
    return iter_->first;
  }

  ///
  inline uint global_index() const
  {
    return distdata_.get_global(iter_->first);
  }

  ///
  inline bool end() const
  {
    return iter_ == distdata_.shared_.end();
  }

  ///
  inline _set<uint> const& adj() const
  {
    return iter_->second;
  }

private:

  DistributedData const& distdata_;
  _map<uint, _set<uint> >::const_iterator iter_;

};

/**
 *  @class  GhostIterator
 *
 *  @brief  Implements an iterator on ghost entities.
 */

class GhostIterator
{

public:

  ///
  GhostIterator(DistributedData const& distdata) :
      distdata_(distdata),
      iter_(distdata_.ghost_.begin())
  {
  }

  ///
  ~GhostIterator()
  {
  }

  ///
  GhostIterator& operator++()
  {
    ++iter_;
    return *this;
  }

  ///
  inline uint index() const
  {
    return iter_->first;
  }

  ///
  inline uint global_index() const
  {
    return distdata_.get_global(iter_->first);
  }

  ///
  inline uint owner() const
  {
    return iter_->second;
  }

  ///
  inline bool end() const
  {
    return iter_ == distdata_.ghost_.end();
  }

  ///
  inline _set<uint> const& adj() const
  {
    return distdata_.shared_.find(iter_->first)->second;
  }

private:

  DistributedData const& distdata_;
  _map<uint, uint>::const_iterator iter_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_DISTRIBUTED_DATA_H */
