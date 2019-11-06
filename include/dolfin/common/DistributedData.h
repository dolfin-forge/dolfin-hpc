// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_COMMON_DISTRIBUTED_DATA_H
#define __DOLFIN_COMMON_DISTRIBUTED_DATA_H

#include <dolfin/common/Distributed.h>

#include <dolfin/common/Array.h>

namespace dolfin
{

class SharedMapping;

/*
 *  @class  DistributedData
 */

class DistributedData : public Distributed<DistributedData>
{

  friend class SharedIterator;
  friend class GhostIterator;
  friend class OwnedIterator;

public:

  ///
  DistributedData(MPI::Communicator& comm = MPI::DOLFIN_COMM);

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

  ///
  friend void swap( DistributedData& a, DistributedData& b )
  {
    using std::swap;

    swap( static_cast<Distributed<DistributedData>&>(a),
          static_cast<Distributed<DistributedData>&>(b) );

    // Swap flags
    swap(a.valid_numbering, b.valid_numbering);
    swap(a.valid_ownership, b.valid_ownership);
    swap(a.valid_adjacency, b.valid_adjacency);

    // Swap attributes
    swap(a.rank_             , b.rank_);
    swap(a.pe_size_          , b.pe_size_);
    swap(a.range_is_set_     , b.range_is_set_);
    swap(a.offset_           , b.offset_);
    swap(a.range_size_       , b.range_size_);
    swap(a.global_size_      , b.global_size_);
    swap(a.finalized_        , b.finalized_);
    swap(a.global_           , b.global_);
    swap(a.local_            , b.local_);
    swap(a.adjacents_        , b.adjacents_);
    swap(a.shared_           , b.shared_);
    swap(a.ghost_            , b.ghost_);
    swap(a.cache_size_       , b.cache_size_);
    swap(a.cached_numbering_ , b.cached_numbering_);
    swap(a.cached_ownership_ , b.cached_ownership_);
    swap(a.shared_mapping_   , b.shared_mapping_);
  }

  /// Finalize the data: validate and set process range + global size
  void finalize();

  /// Assign using other data given mapping from self to other between entities.
  /// Data is finalized.
  void assign(DistributedData const& other, Array<uint> const& mapping);

  /// Return if the distributed data is empty
  bool empty() const;

  /// Return the storage capacity of mappings
  uint capacity() const;

  /// Return whether the data is finalized
  bool is_finalized() const;

  //--- Data bounds -----------------------------------------------------------

  /// Return the offset of the process range
  uint offset() const;

  /// Return the process range size
  uint range_size() const;

  /// Return if the process range is set
  bool range_is_set() const;

  /// Return if the global index is in the process range
  bool in_range(uint global_index) const;

  /// Return if the global index is off the process range
  bool off_range(uint global_index) const;

  /// Return the local data size
  uint local_size() const;

  /// Return the global data size i.e the global sum of number of owned entities
  uint global_size() const;

  /// Set the process range and the global size: if the second is not provided
  /// it is computed by summing the number of owned entities.
  /// Setting range is only possible to a non-finalized distributed data.
  /// If the data is not empty then provided arguments are checked to be
  /// consistent.
  void set_range(uint num_owned, uint num_global = 0);

  /// Set the process local and global size: if the second is not provided
  /// it is computed by summing the number of owned entities.
  /// Setting size is only possible to an empty distributed data and triggers
  /// the creation of cached data structure to avoid use of maps.
  void set_size(uint num_local, uint num_global = 0);

  //--- Numbering -------------------------------------------------------------

  /// Return if the index is a local index
  uint has_local(uint local_index) const;

  /// Return the global index given a local index
  uint get_global(uint local_index) const;

  /// Get the global indices given an array of n local indices
  void get_global(uint n, uint const * local_indices, uint * global_indices) const;

  /// Get the local indices given an array in-place
  inline void get_global(uint n, uint * local_indices) const
  {
    get_global(n, local_indices, local_indices);
  }

  /// Return if the index is a global index
  uint has_global(uint global_index) const;

  /// Return the local index given a global index
  uint get_local(uint global_index) const;

  /// Get the local indices given an array of n local indices
  void get_local(uint n, uint const * global_indices, uint * local_indices) const;

  /// Get the local indices given an array in-place
  inline void get_local(uint n, uint * global_indices) const
  {
    get_local(n, global_indices, global_indices);
  }

  /// Set local-to-global mapping
  void set_map(uint local_index, uint global_index, bool allow_remap = false);

  /// Set local-to-global mapping
  void set_map(Array<uint> const& mapping);

  /// Re-map numbering with given mapping for new local entities numbering
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

  /// Re-map ownership with given mapping for process ranks
  void remap_ownership(Array<uint> const& mapping);

  //--- Shared ---

  /// Return the adjacent set of a shared entity
  _set<uint> const& get_shared_adj(uint local_index) const;

  /// Return the adjacent set of a shared entity
  _set<uint> const* ptr_shared_adj(uint local_index) const;

  /// Return the common adjacent set to an array of shared entities
  void get_common_adj(uint n, uint const indices[], _set<uint>& adjs) const;

  /// Set the entity as shared, the adjacent set is not created.
  /// If the entity is ghost then it stays that way
  void set_shared(uint local_index);

  /// Add a rank as adjacent, this cannot be self
  void set_shared_adj(uint local_index, uint adj);

  /// Set the adjacent set for the given shared entity, this cannot contain self
  void setall_shared_adj(uint local_index, _set<uint> const& adjs);

  /// Return shared entities mapping, only on finalized data
  SharedMapping const& shared_mapping() const;

  /// Re-map the adjacent set
  void remap_shared_adj();

  /// Check shared entities consistency
  void check_shared();

  //--- Ghosts ---

  /// Set the given entity as ghost
  void set_ghost(uint local_index, uint owner);

  /// Check ghost entities consistency
  void check_ghost();

  //---------------------------------------------------------------------------

  //
  void disp() const;

public:

  ///
  bool valid_numbering;
  bool valid_ownership;
  bool valid_adjacency;

private:

  /// Clear all data
  void clear();

  uint rank_;
  uint pe_size_;

  //
  bool range_is_set_;
  uint offset_;
  uint range_size_;

  //
  uint global_size_;

  //
  bool finalized_;

  ///
  typedef _map<uint, uint> IndexMapping;
  IndexMapping global_;
  IndexMapping local_;

  //
  _set<uint> adjacents_;
  typedef _map<uint, _set<uint> > SharedSet;
  SharedSet shared_;
  typedef _map<uint, uint> GhostSet;
  GhostSet ghost_;

  ///
  uint cache_size_;
  uint * cached_numbering_;
  uint * cached_ownership_;

  /// Mapping created on-demand
  mutable SharedMapping * shared_mapping_;

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
  inline uint owner() const
  {
    return distdata_.get_owner(iter_->first);
  }

  ///
  inline bool is_owned() const
  {
    return distdata_.is_owned(iter_->first);
  }

  ///
  inline bool valid() const
  {
    return iter_ != distdata_.shared_.end();
  }

  ///
  inline _set<uint> const& adj() const
  {
    return iter_->second;
  }

  ///
  template <class T>
  inline void adj_enqueue(Array<T> container[], T value) const
  {
    _set<uint> const& a = iter_->second;
    for (_set<uint>::const_iterator it = a.begin(); it != a.end(); ++it)
    {
      container[*it].push_back(value);
    }
  }

private:

  DistributedData const& distdata_;
  DistributedData::SharedSet::const_iterator iter_;

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
  inline bool valid() const
  {
    return iter_ != distdata_.ghost_.end();
  }

  ///
  inline _set<uint> const& adj() const
  {
    return distdata_.shared_.find(iter_->first)->second;
  }

  ///
  template <class T>
  inline void adj_enqueue(Array<T> container[], T value) const
  {
    _set<uint> const& a = distdata_.shared_.find(iter_->first)->second;
    for (_set<uint>::const_iterator it = a.begin(); it != a.end(); ++it)
    {
      container[*it].push_back(value);
    }
  }

private:

  DistributedData const& distdata_;
  DistributedData::GhostSet::const_iterator iter_;

};



/**
 *  @class  OwnedIterator
 *
 *  @brief  Implements an iterator on owned entities for finalized distributed
 *          data only.
 */

class OwnedIterator
{

public:

  ///
  OwnedIterator(DistributedData const& distdata) :
      distdata_(distdata),
      owner_(distdata.cached_ownership_),
      begin_(distdata.cached_numbering_),
      end_(distdata.cached_numbering_ + distdata_.local_size()),
      iter_(begin_)
  {
    if(!distdata.is_finalized())
    {
      error("OwnedIterator : distributed data is not finalized");
    }
  }

  ///
  ~OwnedIterator()
  {
  }

  ///
  OwnedIterator& operator++()
  {
    if (iter_ == end_)
    {
      return *this;
    }
    ++iter_;
    while ((iter_ < end_)
           && (owner_[iter_ - begin_] != distdata_.pe_size_)
           && (owner_[iter_ - begin_] != distdata_.rank_))
    {
      ++iter_;
    }
    return *this;
  }

  ///
  inline uint index() const
  {
    return iter_ - begin_;
  }

  ///
  inline uint global_index() const
  {
    return *iter_;
  }

  ///
  inline uint is_shared() const
  {
    return (owner_[iter_ - begin_] == distdata_.rank_);
  }

  ///
  inline bool end() const
  {
    return iter_ == end_;
  }

private:

  DistributedData const& distdata_;
  uint * const owner_;
  uint * const begin_;
  uint * const end_;
  uint * iter_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_COMMON_DISTRIBUTED_DATA_H */
