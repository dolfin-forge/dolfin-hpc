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
  ~DistributedData() override;

  ///
  auto operator=( DistributedData const & other ) -> DistributedData &;

  ///
  friend void swap( DistributedData & a, DistributedData & b );

  ///
  auto operator==( DistributedData const & other ) const -> bool;

  ///
  auto operator!=( DistributedData const & other ) const -> bool;

  /// Finalize the data: validate and set process range + global size
  void finalize();

  /// Assign using other data given mapping from self to other between entities.
  /// Data is finalized.
  void assign( DistributedData const & other, Array< uint > const & mapping );

  /// Return if the distributed data is empty
  auto empty() const -> bool;

  /// Return the storage capacity of mappings
  auto capacity() const -> uint;

  /// Return whether the data is finalized
  auto is_finalized() const -> bool;

  //--- Data bounds -----------------------------------------------------------

  /// Return the offset of the process range
  auto offset() const -> uint;

  /// Return the process range size
  auto range_size() const -> uint;

  /// Return if the process range is set
  auto range_is_set() const -> bool;

  /// Return if the global index is in the process range
  auto in_range( uint global_index ) const -> bool;

  /// Return if the global index is off the process range
  auto off_range( uint global_index ) const -> bool;

  /// Return the local data size
  auto local_size() const -> uint;

  /// Return the global data size i.e the global sum of number of owned entities
  auto global_size() const -> uint;

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
  auto has_local( uint local_index ) const -> uint;

  /// Return the global index given a local index
  auto get_global( uint local_index ) const -> uint;

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
  auto has_global( uint global_index ) const -> uint;

  /// Return the local index given a global index
  auto get_local( uint global_index ) const -> uint;

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
  auto has_adj_rank( uint rank ) const -> bool;

  /// Return the number of adjacent ranks
  auto num_adj_ranks() const -> uint;

  /// Return the set of adjacent ranks
  auto get_adj_ranks() const -> _set< uint > const &;

  //--- Ownership -------------------------------------------------------------

  /// Return the owner of the entity: self if owned and not self otherwise
  auto get_owner( uint local_index ) const -> uint;

  /// Return if the entity is owned
  auto is_owned( uint local_index ) const -> bool;

  /// Return if the entity is shared: it can be owned or not
  auto is_shared( uint local_index ) const -> bool;

  /// Return if the entity is ghost:  it is shared and not owned
  auto is_ghost( uint local_index ) const -> bool;

  /// Return the number of owned entities
  auto num_owned() const -> uint;

  /// Return the number of shared entities
  auto num_shared() const -> uint;

  /// Return the number of ghost entities
  auto num_ghost() const -> uint;

  /// Re-map ownership with given mapping for process ranks
  void remap_ownership( Array< uint > const & mapping );

  //--- Shared ---

  /// Return the adjacent set of a shared entity
  auto get_shared_adj( uint local_index ) const -> _set< uint > const &;

  /// Return the adjacent set of a shared entity
  auto ptr_shared_adj( uint local_index ) const -> _set< uint > const *;

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
  auto shared_mapping() const -> SharedMapping const &;

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
  bool valid_numbering{ false };

private:
  /// Clear all data
  void clear();

  uint rank_;
  uint pe_size_;

  //
  bool range_is_set_{ false };
  uint offset_{ 0 };
  uint range_size_{ 0 };

  //
  uint global_size_{ 0 };

  //
  bool finalized_{ false };

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
  mutable SharedMapping * shared_mapping_{};
};
//-----------------------------------------------------------------------------
inline auto DistributedData::operator==( DistributedData const & other ) const -> bool
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
inline auto DistributedData::operator!=( DistributedData const & other ) const -> bool
{
  return !( *this == other );
}
//-----------------------------------------------------------------------------
inline auto DistributedData::empty() const -> bool
{
  return ( local_.empty() and global_.empty() and shared_.empty()
           and ghost_.empty() );
}
//-----------------------------------------------------------------------------
inline auto DistributedData::capacity() const -> uint
{
  return local_.size();
}
//-----------------------------------------------------------------------------
inline auto DistributedData::is_finalized() const -> bool
{
  return finalized_;
}
//-----------------------------------------------------------------------------
inline auto DistributedData::offset() const -> uint
{
  return offset_;
}
//-----------------------------------------------------------------------------
inline auto DistributedData::range_size() const -> uint
{
  return range_size_;
}
//-----------------------------------------------------------------------------
inline auto DistributedData::range_is_set() const -> bool
{
  return range_is_set_;
}
//-----------------------------------------------------------------------------
inline auto DistributedData::in_range( uint global_index ) const -> bool
{
  dolfin_assert( global_size_ > 0 );
  dolfin_assert( offset_ + range_size_ <= global_size_ );
  return ( offset_ <= global_index && global_index < offset_ + range_size_ );
}
//-----------------------------------------------------------------------------
inline auto DistributedData::off_range( uint global_index ) const -> bool
{
  dolfin_assert( global_size_ > 0 );
  dolfin_assert( offset_ + range_size_ <= global_size_ );
  return ( global_index < offset_ || offset_ + range_size_ <= global_index );
}
//-----------------------------------------------------------------------------
inline auto DistributedData::local_size() const -> uint
{
  // If local size is not known, return current size, otherwise return
  return ( cached_numbering_.empty() ? local_.size()
                                     : cached_numbering_.size() );
}
//-----------------------------------------------------------------------------
inline auto DistributedData::global_size() const -> uint
{
  if ( global_size_ < local_.size() )
  {
    error( "DistributedData : global size has not been set or is invalid" );
  }
  return global_size_;
}
//-----------------------------------------------------------------------------
inline auto DistributedData::has_local( uint local_index ) const -> uint
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
inline auto DistributedData::get_global( uint local_index ) const -> uint
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
inline auto DistributedData::has_global( uint global_index ) const -> uint
{
  return ( local_.count( global_index ) > 0 );
}
//-----------------------------------------------------------------------------
inline auto DistributedData::get_local( uint global_index ) const -> uint
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
inline auto DistributedData::has_adj_rank( uint rank ) const -> bool
{
  return ( adjacents_.count( rank ) > 0 );
}
//-----------------------------------------------------------------------------
inline auto DistributedData::num_adj_ranks() const -> uint
{
  return adjacents_.size();
}
//-----------------------------------------------------------------------------
inline auto DistributedData::get_adj_ranks() const -> _set< uint > const &
{
  return adjacents_;
}
//-----------------------------------------------------------------------------
inline auto DistributedData::get_owner( uint local_index ) const -> uint
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
inline auto DistributedData::is_owned( uint local_index ) const -> bool
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
inline auto DistributedData::is_shared( uint local_index ) const -> bool
{
  if ( not cached_ownership_.empty() )
  {
    dolfin_assert( local_index < cached_ownership_.size() );
    return ( cached_ownership_[local_index] < pe_size_ );
  }
  return ( shared_.count( local_index ) > 0 );
}
//-----------------------------------------------------------------------------
inline auto DistributedData::is_ghost( uint local_index ) const -> bool
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
inline auto DistributedData::num_owned() const -> uint
{
  dolfin_assert( local_.size() >= ghost_.size() );
  return ( local_.size() - ghost_.size() );
}
//-----------------------------------------------------------------------------
inline auto DistributedData::num_shared() const -> uint
{
  return ( shared_.size() );
}
//-----------------------------------------------------------------------------
inline auto DistributedData::num_ghost() const -> uint
{
  return ( ghost_.size() );
}
//-----------------------------------------------------------------------------
inline auto
  DistributedData::get_shared_adj( uint local_index ) const -> _set< uint > const &
{
  dolfin_assert( shared_.count( local_index ) > 0 );
  return shared_.find( local_index )->second;
}
//-----------------------------------------------------------------------------
inline auto
  DistributedData::ptr_shared_adj( uint local_index ) const -> _set< uint > const *
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
