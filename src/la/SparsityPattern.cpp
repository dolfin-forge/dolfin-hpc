// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/la/SparsityPattern.h>

#include <dolfin/common/timing.h>
#include <dolfin/log/LogStream.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

#include <iostream>

namespace dolfin
{

//-----------------------------------------------------------------------------
SparsityPattern::SparsityPattern( size_t         rank,
                                  size_t const * dim,
                                  size_t const * range /* = nullptr */ )
  : rank_( 0 )
  , dim_( nullptr )
  , range_( nullptr )
  , local_range_( nullptr )
  , initialized_( false )
  , finalized_( false )
  , blocked_( false )
  , distributed_( false )
  , d_entries_( nullptr )
  , d_count_( 0 )
  , o_entries_( nullptr )
  , o_count_( 0 )
{
  init( rank, dim, range );
}
//-----------------------------------------------------------------------------
SparsityPattern::~SparsityPattern()
{
  clear();
}
//-----------------------------------------------------------------------------
void SparsityPattern::clear()
{
  r_entries_.clear();
  o_count_ = 0;
  delete[] o_entries_;
  o_entries_ = nullptr;
  d_count_   = 0;
  delete[] d_entries_;
  d_entries_   = nullptr;
  initialized_ = false;
  finalized_   = false;
  blocked_     = false;
  distributed_ = false;
  for ( size_t i = 0; i < rank_; ++i )
  {
    delete[] range_[i];
  }
  delete[] range_;
  range_ = nullptr;
  delete[] local_range_;
  local_range_ = nullptr;
  delete[] dim_;
  dim_  = nullptr;
  rank_ = 0;
}
//-----------------------------------------------------------------------------
void SparsityPattern::init( size_t         rank,
                            size_t const * dim,
                            size_t const * range /* = nullptr */ )
{
  if ( initialized_ )
  {
    error( "SparsityPattern : pattern has already been initialized" );
  }
  if ( rank == 0 )
  {
    error( "SparsityPattern : unimplemented for rank = 0" );
  }
  if ( rank > 2 )
  {
    error( "SparsityPattern : unimplemented for rank > 2" );
  }

  // Set tensor rank and global dimensions
  rank_ = rank;
  dim_  = new size_t[rank];
  std::copy( dim, dim + rank, dim_ );

  // Set the sparsity pattern as distributed only if ranges are provided and not
  // all equal to the global dimensions
  distributed_ =
    ( range != nullptr ) && ( !std::equal( dim, dim + rank, range ) );
  if ( distributed_ )
  {
    size_t pe_size = dolfin::MPI::size();
    size_t pe_rank = dolfin::MPI::rank();
    range_         = new size_t *[rank];
    local_range_   = new size_t *[rank];
    for ( size_t i = 0; i < rank; ++i )
    {
      range_[i]       = new size_t[pe_size + 1];
      local_range_[i] = &range_[i][pe_rank];
      // Collect all the range size in the array then sum to compute offsets
      // Previous code created a temporary array, why not do it in-place.
      MPI::all_gather(
        ( size_t * ) &range[i], 1, ( size_t * ) &range_[i][1], 1 );
      range_[i][0] = 0;
      for ( size_t j = 0; j < pe_size; ++j )
      {
        range_[i][j + 1] += range_[i][j];
      }
    }
  }
  else
  {
    // Since the pattern is not distributed only store the process range, this
    // means that calls to functions assuming distributed pattern will be
    // invalid.
    range_       = new size_t *[rank];
    local_range_ = new size_t *[rank];
    for ( size_t i = 0; i < rank; ++i )
    {
      range_[i]       = new size_t[2];
      local_range_[i] = &range_[i][0];
      range_[i][0]    = 0;
      range_[i][1]    = dim[i];
    }
  }

  // This data structure contains set of non-zero column for each row
  // in the process range: since this range has a given size and every
  // row has at least a non-zero entry
  d_entries_ = new _ordered_set< size_t >[this->size( 0 )];
  if ( distributed_ )
  {
    o_entries_ = new _ordered_set< size_t >[this->size( 0 )];
  }

  //
  initialized_ = true;
}
//-----------------------------------------------------------------------------
void SparsityPattern::insert( size_t const * num, size_t const * const * idx )
{
  for ( size_t i = 0; i < num[0]; ++i )
  {
    // Row is in the process range
    if ( local_range_[0][0] <= idx[0][i] && idx[0][i] < local_range_[0][1] )
    {
      for ( size_t j = 0; j < num[1]; ++j )
      {
        // Column is in the process range: diagonal portion
        if ( local_range_[1][0] <= idx[1][j] && idx[1][j] < local_range_[1][1] )
        {
          d_entries_[idx[0][i] - local_range_[0][0]].insert( idx[1][j] );
          ++d_count_;
        }
        // Off-diagonal portion
        else
        {
          dolfin_assert( distributed_ );
          o_entries_[idx[0][i] - local_range_[0][0]].insert( idx[1][j] );
          ++o_count_;
        }
      }
    }
    else
    {
      dolfin_assert( distributed_ );
      for ( size_t j = 0; j < num[1]; ++j )
      {
        r_entries_[idx[0][i]].insert( idx[1][j] );
      }
    }
  }
}
//-----------------------------------------------------------------------------
void SparsityPattern::numNonZeroPerRow( size_t nzrow[] ) const
{
  if ( rank_ != 2 )
  {
    error(
      "SparsityPattern : non-zero entries per row can only be computed for "
      "matrices." );
  }
  if ( d_count_ == 0 && o_count_ == 0 )
  {
    error( "SparsityPattern : pattern has not been computed." );
  }

  size_t const num_rows = this->size( 0 );
  for ( size_t i = 0; i < num_rows; ++i )
  {
    nzrow[i] = d_entries_[i].size();
  }

  if ( o_entries_ != nullptr )
  {
    for ( size_t i = 0; i < num_rows; ++i )
    {
      nzrow[i] += o_entries_[i].size();
    }
  }
}
//-----------------------------------------------------------------------------
void SparsityPattern::numNonZeroPerRow( size_t p_rank,
                                        size_t d_nzrow[],
                                        size_t o_nzrow[] ) const
{
  if ( rank_ != 2 )
  {
    error( "Non-zero entries per row can be computed for matrices only." );
  }
  if ( d_count_ == 0 && o_count_ == 0 )
  {
    error( "SparsityPattern : pattern has not been computed." );
  }

  /// All information is in the local range
  if ( p_rank == MPI::rank() )
  {
    size_t const num_rows = this->size( 0 );
    for ( size_t i = 0; i < num_rows; ++i )
    {
      d_nzrow[i] = d_entries_[i].size();
      o_nzrow[i] = o_entries_[i].size();
    }
  }
  else
  {
    size_t const r0 = range_[0][p_rank];
    size_t const r1 = range_[0][p_rank + 1];
    size_t const c0 = range_[1][p_rank];
    size_t const c1 = range_[1][p_rank + 1];
    std::fill_n( d_nzrow, r1 - r0, 0 );
    std::fill_n( o_nzrow, r1 - r0, 0 );
    for ( _ordered_map< size_t, _ordered_set< size_t > >::const_iterator it =
            r_entries_.find( r0 );
          it->first < r1;
          ++it )
    {
      for ( _ordered_set< size_t >::const_iterator c = it->second.begin();
            c != it->second.end();
            ++c )
      {
        if ( ( c0 <= *c ) && ( *c < c1 ) )
        {
          ++d_nzrow[it->first - r0];
        }
        else
        {
          ++o_nzrow[it->first - r0];
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
auto SparsityPattern::numNonZero() const -> size_t
{
  if ( rank_ != 2 )
  {
    error( "SparsityPattern : total non-zeros entries can only be computed for "
           "matrices." );
  }
  if ( d_count_ == 0 && o_count_ == 0 )
  {
    error( "SparsityPattern : pattern has not been computed." );
  }

  // Compute total number of nonzeros per row
  size_t       nz       = 0;
  size_t const num_rows = this->size( 0 );
  for ( size_t i = 0; i < num_rows; ++i )
  {
    nz += d_entries_[i].size();
  }

  return nz;
}
//-----------------------------------------------------------------------------
void SparsityPattern::disp() const
{
  section( "SparsityPattern" );
  message( "distributed : %u", distributed_ );
  message( "rank        : %u", rank_ );
  begin( "dim         :" );
  for ( size_t i = 0; i < rank_; ++i )
  {
    message( "%u : %u", i, dim_[i] );
  }
  end();
  message( "range       : [ %8u, %8u [", local_range_[0], local_range_[1] );
  end();
}
//-----------------------------------------------------------------------------
void SparsityPattern::apply()
{
  finalized_ = true;

  if ( !distributed_ )
    return;

#ifdef DOLFIN_HAVE_MPI

  message(
    1, "SparsityPattern : apply with %u remote rows", r_entries_.size() );
  tic();

  /// Collect entries per owner
  size_t const rank    = MPI::rank();
  size_t const pe_size = MPI::size();
  size_t       owner   = 0;
  size_t       sendmax = 0;

  std::vector< std::vector< size_t > > sendbuf( pe_size );

  for ( _ordered_map< size_t, _ordered_set< size_t > >::const_iterator it =
          r_entries_.begin(); it != r_entries_.end(); ++it )
  {
    // Increment owner when jumping to another range
    while ( range_[0][owner + 1] <= it->first )
    {
      sendmax = std::max( sendmax, sendbuf[owner].size() );
      ++owner;
    }

    dolfin_assert( range_[0][owner] <= it->first );
    dolfin_assert( it->first < range_[0][owner + 1] );

    // Data packet [ global index, number of entries, [ indices ] ]
    sendbuf[owner].push_back( it->first );
    sendbuf[owner].push_back( it->second.size() );
    sendbuf[owner].insert( sendbuf[owner].end(),
                           it->second.begin(),
                           it->second.end() );
  }

  size_t recvmax = std::max( sendmax, sendbuf[owner].size() );
  MPI::all_reduce_in_place< MPI::max >( recvmax );
  if ( recvmax == 0 )
  {
    return;
  }

  ///
  size_t const r0 = local_range_[0][0];
  size_t const r1 = local_range_[0][1];
  size_t const c0 = local_range_[1][0];
  size_t const c1 = local_range_[1][1];

  /// Exchange entries and add to diagonal and off-diagonal data structures
  std::vector< size_t > recvbuf( recvmax );
  for ( size_t j = 1; j < pe_size; ++j )
  {
    size_t src = ( rank - j + pe_size ) % pe_size;
    size_t dst = ( rank + j ) % pe_size;

    int recv_count = MPI::sendrecv( sendbuf[dst], dst, recvbuf, src, 1 );

    for ( int k = 0; k < recv_count; )
    {
      size_t const irow = recvbuf[k++];
      if ( irow < r0 || r1 <= irow )
      {
        error( "SparsityPattern : invalid global row number received %u",
               irow );
      }
      size_t const ncol = recvbuf[k++];
      for ( size_t p = 0; p < ncol; ++p, ++k )
      {
        if ( c0 <= recvbuf[k] && recvbuf[k] < c1 )
        {
          d_entries_[irow - r0].insert( recvbuf[k] );
          ++d_count_;
        }
        else
        {
          o_entries_[irow - r0].insert( recvbuf[k] );
          ++o_count_;
        }
      }
    }
  }

  tocd( 1 );

#endif
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
