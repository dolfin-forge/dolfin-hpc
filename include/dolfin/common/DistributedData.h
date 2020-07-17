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

class DistributedData : public Distributed< DistributedData >
{

  friend class GhostIterator;
  friend class OwnedIterator;
  friend class SharedIterator;

public:
  using IndexMapping = _map< uint, uint >;
  using SharedSet    = _map< uint, _set< uint > >;
  using GhostSet     = _map< uint, uint >;

public:
  ///
  DistributedData( MPI::Communicator & comm = MPI::DOLFIN_COMM );

  ///
  DistributedData( DistributedData const & other );

  ///
  ~DistributedData();

  ///
  DistributedData & operator=( DistributedData const & other );

  ///
  friend void swap( DistributedData & a, DistributedData & b );

  ///
  bool operator==( DistributedData const & other ) const;

  ///
  bool operator!=( DistributedData const & other ) const;

  /// Finalize the data: validate and set process range + global size
  void finalize();

  /// Assign using other data given mapping from self to other between entities.
  /// Data is finalized.
  void assign( DistributedData const & other, Array< uint > const & mapping );

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
  bool in_range( uint global_index ) const;

  /// Return if the global index is off the process range
  bool off_range( uint global_index ) const;

  /// Return the local data size
  uint local_size() const;

  /// Return the global data size i.e the global sum of number of owned entities
  uint global_size() const;

  /// Set the process range and the global size: if the second is not provided
  /// it is computed by summing the number of owned entities.
  /// Setting range is only possible to a non-finalized distributed data.
  /// If the data is not empty then provided arguments are checked to be
  /// consistent.
  void set_range( uint num_owned, uint num_global = 0 );

  /// Set the process local and global size: if the second is not provided
  /// it is computed by summing the number of owned entities.
  /// Setting size is only possible to an empty distributed data and triggers
  /// the creation of cached data structure to avoid use of maps.
  void set_size( uint num_local, uint num_global = 0 );

  //--- Numbering -------------------------------------------------------------

  /// Return if the index is a local index
  uint has_local( uint local_index ) const;

  /// Return the global index given a local index
  uint get_global( uint local_index ) const;

  /// Get the global indices given an array of n local indices
  void get_global( uint         n,
                   uint const * local_indices,
                   uint *       global_indices ) const;

  /// Get the local indices given an array in-place
  inline void get_global( uint n, uint * local_indices ) const
  {
    get_global( n, local_indices, local_indices );
  }

  /// Return if the index is a global index
  uint has_global( uint global_index ) const;

  /// Return the local index given a global index
  uint get_local( uint global_index ) const;

  /// Get the local indices given an array of n local indices
  void get_local( uint         n,
                  uint const * global_indices,
                  uint *       local_indices ) const;

  /// Get the local indices given an array in-place
  inline void get_local( uint n, uint * global_indices ) const
  {
    get_local( n, global_indices, global_indices );
  }

  /// Set local-to-global mapping
  void set_map( uint local_index, uint global_index, bool allow_remap = false );

  /// Set local-to-global mapping
  void set_map( Array< uint > const & mapping );

  /// Re-map numbering with given mapping for new local entities numbering
  void remap_numbering( Array< uint > const & mapping );

  /// Re-index global indices to have contiguous numbering of owned entities
  void renumber_global();

  //--- Adjacency -------------------------------------------------------------

  /// Return whether the given rank is adjacent
  bool has_adj_rank( uint rank ) const;

  /// Return the number of adjacent ranks
  uint num_adj_ranks() const;

  /// Return the set of adjacent ranks
  _set< uint > const & get_adj_ranks() const;

  //--- Ownership -------------------------------------------------------------

  /// Return the owner of the entity: self if owned and not self otherwise
  uint get_owner( uint local_index ) const;

  /// Return if the entity is owned
  bool is_owned( uint local_index ) const;

  /// Return if the entity is shared: it can be owned or not
  bool is_shared( uint local_index ) const;

  /// Return if the entity is ghost:  it is shared and not owned
  bool is_ghost( uint local_index ) const;

  /// Return the number of owned entities
  uint num_owned() const;

  /// Return the number of shared entities
  uint num_shared() const;

  /// Return the number of ghost entities
  uint num_ghost() const;

  /// Re-map ownership with given mapping for process ranks
  void remap_ownership( Array< uint > const & mapping );

  //--- Shared ---

  /// Return the adjacent set of a shared entity
  _set< uint > const & get_shared_adj( uint local_index ) const;

  /// Return the adjacent set of a shared entity
  _set< uint > const * ptr_shared_adj( uint local_index ) const;

  /// Return the common adjacent set to an array of shared entities
  void get_common_adj( uint n, Array< uint > const & indices,
                       _set< uint > & adjs ) const;

  /// Set the entity as shared, the adjacent set is not created.
  /// If the entity is ghost then it stays that way
  void set_shared( uint local_index );

  /// Add a rank as adjacent, this cannot be self
  void set_shared_adj( uint local_index, uint adj );

  /// Set the adjacent set for the given shared entity, this cannot contain self
  void setall_shared_adj( uint local_index, _set< uint > const & adjs );

  /// Return shared entities mapping, only on finalized data
  SharedMapping const & shared_mapping() const;

  /// Re-map the adjacent set
  void remap_shared_adj();

  /// Check shared entities consistency
  void check_shared();

  //--- Ghosts ---

  /// Set the given entity as ghost
  void set_ghost( uint local_index, uint owner );

  /// Check ghost entities consistency
  void check_ghost();

  //---------------------------------------------------------------------------

  //
  void disp() const;

public:
  ///
  bool valid_numbering;

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
  IndexMapping global_;
  IndexMapping local_;

  //
  _set< uint > adjacents_;
  SharedSet    shared_;
  GhostSet     ghost_;

  ///
  Array< uint > cached_numbering_;
  Array< uint > cached_ownership_;

  /// Mapping created on-demand
  mutable SharedMapping * shared_mapping_;
};
//-----------------------------------------------------------------------------
inline bool DistributedData::operator==( DistributedData const & other ) const
{
  if ( rank_ != other.rank_ )
    return false;

  if ( pe_size_ != other.pe_size_ )
    return false;

  if ( range_is_set_ != other.range_is_set_ )
    return false;

  if ( offset_ != other.offset_ )
    return false;

  if ( range_size_ != other.range_size_ )
    return false;

  if ( global_size_ != other.global_size_ )
    return false;

  if ( finalized_ != other.finalized_ )
    return false;

  if ( cached_numbering_ != other.cached_numbering_ )
    return false;

  if ( cached_ownership_ != other.cached_ownership_ )
    return false;

  return true;
}
//-----------------------------------------------------------------------------
inline bool DistributedData::operator!=( DistributedData const & other ) const
{
  return !( *this == other );
}
//-----------------------------------------------------------------------------
inline bool DistributedData::empty() const
{
  return ( local_.empty() and global_.empty() and shared_.empty()
           and ghost_.empty() );
}
//-----------------------------------------------------------------------------
inline uint DistributedData::capacity() const
{
  return local_.size();
}
//-----------------------------------------------------------------------------
inline bool DistributedData::is_finalized() const
{
  return finalized_;
}
//-----------------------------------------------------------------------------
inline uint DistributedData::offset() const
{
  return offset_;
}
//-----------------------------------------------------------------------------
inline uint DistributedData::range_size() const
{
  return range_size_;
}
//-----------------------------------------------------------------------------
inline bool DistributedData::range_is_set() const
{
  return range_is_set_;
}
//-----------------------------------------------------------------------------
inline bool DistributedData::in_range( uint global_index ) const
{
  dolfin_assert( global_size_ > 0 );
  dolfin_assert( offset_ + range_size_ <= global_size_ );
  return ( offset_ <= global_index && global_index < offset_ + range_size_ );
}
//-----------------------------------------------------------------------------
inline bool DistributedData::off_range( uint global_index ) const
{
  dolfin_assert( global_size_ > 0 );
  dolfin_assert( offset_ + range_size_ <= global_size_ );
  return ( global_index < offset_ || offset_ + range_size_ <= global_index );
}
//-----------------------------------------------------------------------------
inline uint DistributedData::local_size() const
{
  // If local size is not known, return current size, otherwise return
  return ( cached_numbering_.empty() ? local_.size()
                                     : cached_numbering_.size() );
}
//-----------------------------------------------------------------------------
inline uint DistributedData::global_size() const
{
  if ( global_size_ < local_.size() )
  {
    error( "DistributedData : global size has not been set or is invalid" );
  }
  return global_size_;
}
//-----------------------------------------------------------------------------
inline uint DistributedData::has_local( uint local_index ) const
{
  if ( not cached_numbering_.empty() )
  {
    dolfin_assert( global_.size() == 0 );
    dolfin_assert( local_index < cached_numbering_.size() );
    return cached_numbering_[local_index] != DOLFIN_UINT_UNDEF;
  }
  return ( global_.count( local_index ) > 0 );
}
//-----------------------------------------------------------------------------
inline uint DistributedData::get_global( uint local_index ) const
{
  if ( not cached_numbering_.empty() )
  {
    dolfin_assert( global_.size() == 0 );
    dolfin_assert( local_index < cached_numbering_.size() );
    dolfin_assert( cached_numbering_[local_index] != DOLFIN_UINT_UNDEF );
    return cached_numbering_[local_index];
  }
  dolfin_assert( global_.count( local_index ) > 0 );
  return global_.find( local_index )->second;
}
//-----------------------------------------------------------------------------
inline void DistributedData::get_global( uint         n,
                                         uint const * local_indices,
                                         uint *       global_indices ) const
{
  if ( not cached_numbering_.empty() )
  {
    dolfin_assert( global_.size() == 0 );
    for ( uint i = 0; i < n; ++i )
    {
      dolfin_assert( local_indices[i] < cached_numbering_.size() );
      dolfin_assert( cached_numbering_[local_indices[i]] != DOLFIN_UINT_UNDEF );
      global_indices[i] = cached_numbering_[local_indices[i]];
    }
  }
  else
  {
    for ( uint i = 0; i < n; ++i )
    {
      dolfin_assert( global_.count( local_indices[i] ) > 0 );
      global_indices[i] = global_.find( local_indices[i] )->second;
    }
  }
}
//-----------------------------------------------------------------------------
inline uint DistributedData::has_global( uint global_index ) const
{
  return ( local_.count( global_index ) > 0 );
}
//-----------------------------------------------------------------------------
inline uint DistributedData::get_local( uint global_index ) const
{
  dolfin_assert( local_.count( global_index ) > 0 );
  return local_.find( global_index )->second;
}
//-----------------------------------------------------------------------------
inline void DistributedData::get_local( uint         n,
                                        uint const * global_indices,
                                        uint *       local_indices ) const
{
  for ( uint i = 0; i < n; ++i )
  {
    dolfin_assert( local_.count( global_indices[i] ) > 0 );
    local_indices[i] = local_.find( global_indices[i] )->second;
  }
}
//-----------------------------------------------------------------------------
inline bool DistributedData::has_adj_rank( uint rank ) const
{
  return ( adjacents_.count( rank ) > 0 );
}
//-----------------------------------------------------------------------------
inline uint DistributedData::num_adj_ranks() const
{
  return adjacents_.size();
}
//-----------------------------------------------------------------------------
inline _set< uint > const & DistributedData::get_adj_ranks() const
{
  return adjacents_;
}
//-----------------------------------------------------------------------------
inline uint DistributedData::get_owner( uint local_index ) const
{
  if ( not cached_ownership_.empty() )
  {
    dolfin_assert( local_index < cached_ownership_.size() );
    return ( cached_ownership_[local_index] == pe_size_
               ? rank_
               : cached_ownership_[local_index] );
  }
  GhostSet::const_iterator it = ghost_.find( local_index );
  if ( it == ghost_.end() )
  {
    return rank_;
  }
  return it->second;
}
//-----------------------------------------------------------------------------
inline bool DistributedData::is_owned( uint local_index ) const
{
  if ( not cached_ownership_.empty() )
  {
    dolfin_assert( local_index < cached_ownership_.size() );
    return ( cached_ownership_[local_index] == pe_size_
             || cached_ownership_[local_index] == rank_ );
  }
  return ( ghost_.count( local_index ) == 0 );
}
//-----------------------------------------------------------------------------
inline bool DistributedData::is_shared( uint local_index ) const
{
  if ( not cached_ownership_.empty() )
  {
    dolfin_assert( local_index < cached_ownership_.size() );
    return ( cached_ownership_[local_index] < pe_size_ );
  }
  return ( shared_.count( local_index ) > 0 );
}
//-----------------------------------------------------------------------------
inline bool DistributedData::is_ghost( uint local_index ) const
{
  if ( not cached_ownership_.empty() )
  {
    dolfin_assert( local_index < cached_ownership_.size() );
    return ( cached_ownership_[local_index] < pe_size_
             && cached_ownership_[local_index] != rank_ );
  }
  return ( ghost_.count( local_index ) > 0 );
}
//-----------------------------------------------------------------------------
inline uint DistributedData::num_owned() const
{
  dolfin_assert( local_.size() >= ghost_.size() );
  return ( local_.size() - ghost_.size() );
}
//-----------------------------------------------------------------------------
inline uint DistributedData::num_shared() const
{
  return ( shared_.size() );
}
//-----------------------------------------------------------------------------
inline uint DistributedData::num_ghost() const
{
  return ( ghost_.size() );
}
//-----------------------------------------------------------------------------
inline _set< uint > const &
  DistributedData::get_shared_adj( uint local_index ) const
{
  dolfin_assert( shared_.count( local_index ) > 0 );
  return shared_.find( local_index )->second;
}
//-----------------------------------------------------------------------------
inline _set< uint > const *
  DistributedData::ptr_shared_adj( uint local_index ) const
{
  if ( not cached_ownership_.empty() )
  {
    if ( cached_ownership_[local_index] < pe_size_ )
      return &shared_.find( local_index )->second;
    else
      return nullptr;
  }
  SharedSet::const_iterator it = shared_.find( local_index );
  return ( it == shared_.end() ? nullptr : &it->second );
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_COMMON_DISTRIBUTED_DATA_H */
