// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/la/SparsityPattern.h>

#include <dolfin/common/timing.h>
#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>
#include <dolfin/main/MPI.h>

#include <iostream>

namespace dolfin
{

//-----------------------------------------------------------------------------
SparsityPattern::SparsityPattern(uint rank, uint const * dim,
                                 uint const * range /* = NULL */) :
    rank_(0),
    dim_(NULL),
    range_(NULL),
    local_range_(NULL),
    initialized_(false),
    finalized_(false),
    blocked_(false),
    distributed_(false),
    d_entries_(NULL),
    d_count_(0),
    o_entries_(NULL),
    o_count_(0)
{
  init(rank, dim, range);
}
//-----------------------------------------------------------------------------
SparsityPattern::SparsityPattern() :
    rank_(0),
    dim_(NULL),
    range_(NULL),
    local_range_(NULL),
    initialized_(false),
    finalized_(false),
    blocked_(false),
    distributed_(false),
    d_entries_(NULL),
    d_count_(0),
    o_entries_(NULL),
    o_count_(0)
{
  // Do nothing
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
  o_entries_ = NULL;
  d_count_ = 0;
  delete[] d_entries_;
  d_entries_ = NULL;
  initialized_ = false;
  finalized_ = false;
  blocked_ = false;
  distributed_ = false;
  for (uint i = 0; i < rank_; ++i)
  {
    delete[] range_[i];
  }
  delete[] range_;
  range_ = NULL;
  delete[] local_range_;
  local_range_ = NULL;
  delete[] dim_;
  dim_ = NULL;
  rank_ = 0;
}
//-----------------------------------------------------------------------------
void SparsityPattern::init(uint rank, uint const * dim,
                           uint const * range /* = NULL */)
{
  if (initialized_)
  {
    error("SparsityPattern : pattern has already been initialized");
  }
  if (rank == 0)
  {
    error("SparsityPattern : unimplemented for rank = 0");
  }
  if (rank > 2)
  {
    error("SparsityPattern : unimplemented for rank > 2");
  }

  // Set tensor rank and global dimensions
  rank_ = rank;
  dim_ = new uint[rank];
  std::copy(dim, dim + rank, dim_);

  // Set the sparsity pattern as distributed only if ranges are provided and not
  // all equal to the global dimensions
  distributed_ = (range != NULL) && (!std::equal(dim, dim + rank, range));
  if (distributed_)
  {
    uint pe_size = dolfin::MPI::size();
    uint pe_rank = dolfin::MPI::rank();
    range_ = new uint*[rank];
    local_range_ = new uint*[rank];
    for (uint i = 0; i < rank; ++i)
    {
      range_[i] = new uint[pe_size + 1];
      local_range_[i] = &range_[i][pe_rank];
      // Collect all the range size in the array then sum to compute offsets
      // Previous code created a temporary array, why not do it in-place.
      MPI::all_gather( (uint*) &range[i], 1, (uint*) &range_[i][1], 1 );
      range_[i][0] = 0;
      for (uint j = 0; j < pe_size; ++j)
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
    range_ = new uint*[rank];
    local_range_ = new uint*[rank];
    for (uint i = 0; i < rank; ++i)
    {
      range_[i] = new uint[2];
      local_range_[i] = &range_[i][0];
      range_[i][0] = 0;
      range_[i][1] = dim[i];
    }
  }

  // This data structure contains set of non-zero column for each row
  // in the process range: since this range has a given size and every
  // row has at least a non-zero entry
  d_entries_ = new _ordered_set<uint>[this->size(0)];
  if (distributed_)
  {
    o_entries_ = new _ordered_set<uint>[this->size(0)];
  }

  //
  initialized_ = true;
}
//-----------------------------------------------------------------------------
void SparsityPattern::insert(uint const * num, uint const * const * idx)
{
  for (uint i = 0; i < num[0]; ++i)
  {
    // Row is in the process range
    if (local_range_[0][0] <= idx[0][i] && idx[0][i] < local_range_[0][1])
    {
      for (uint j = 0; j < num[1]; ++j)
      {
        // Column is in the process range: diagonal portion
        if (local_range_[1][0] <= idx[1][j] && idx[1][j] < local_range_[1][1])
        {
          d_entries_[idx[0][i] - local_range_[0][0]].insert(idx[1][j]);
          ++d_count_;
        }
        // Off-diagonal portion
        else
        {
          dolfin_assert(distributed_);
          o_entries_[idx[0][i] - local_range_[0][0]].insert(idx[1][j]);
          ++o_count_;
        }
      }
    }
    else
    {
      dolfin_assert(distributed_);
      for (uint j = 0; j < num[1]; ++j)
      {
        r_entries_[idx[0][i]].insert(idx[1][j]);
      }
    }
  }
}
//-----------------------------------------------------------------------------
void SparsityPattern::numNonZeroPerRow(uint nzrow[]) const
{
  if (rank_ != 2)
  {
    error("SparsityPattern : non-zero entries per row can only be computed for "
          "matrices.");
  }
  if (d_count_ == 0 && o_count_ == 0)
  {
    error("SparsityPattern : pattern has not been computed.");
  }

  uint const num_rows = this->size(0);
  for (uint i = 0; i < num_rows; ++i)
  {
    nzrow[i] = d_entries_[i].size();
  }

  if (o_entries_ != NULL)
  {
    for (uint i = 0; i < num_rows; ++i)
    {
      nzrow[i] += o_entries_[i].size();
    }
  }
}
//-----------------------------------------------------------------------------
void SparsityPattern::numNonZeroPerRow(uint p_rank, uint d_nzrow[],
                                       uint o_nzrow[]) const
{
  if (rank_ != 2)
  {
    error("Non-zero entries per row can be computed for matrices only.");
  }
  if (d_count_ == 0 && o_count_ == 0)
  {
    error("SparsityPattern : pattern has not been computed.");
  }

  /// All information is in the local range
  if (p_rank == MPI::rank())
  {
    uint const num_rows = this->size(0);
    for (uint i = 0; i < num_rows; ++i)
    {
      d_nzrow[i] = d_entries_[i].size();
      o_nzrow[i] = o_entries_[i].size();
    }
  }
  else
  {
    uint const r0 = range_[0][p_rank];
    uint const r1 = range_[0][p_rank + 1];
    uint const c0 = range_[1][p_rank];
    uint const c1 = range_[1][p_rank + 1];
    std::fill_n(d_nzrow, r1 - r0, 0);
    std::fill_n(o_nzrow, r1 - r0, 0);
    for (_ordered_map<uint, _ordered_set<uint> >::const_iterator it =
         r_entries_.find( r0); it->first < r1; ++it)
    {
      for (_ordered_set<uint>::const_iterator c = it->second.begin();
          c != it->second.end(); ++c)
      {
        if ((c0 <= *c) && (*c < c1))
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
uint SparsityPattern::numNonZero() const
{
  if (rank_ != 2)
  {
    error("SparsityPattern : total non-zeros entries can only be computed for "
          "matrices.");
  }
  if (d_count_ == 0 && o_count_ == 0)
  {
    error("SparsityPattern : pattern has not been computed.");
  }

  // Compute total number of nonzeros per row
  uint nz = 0;
  uint const num_rows = this->size(0);
  for (uint i = 0; i < num_rows; ++i)
  {
    nz += d_entries_[i].size();
  }

  return nz;
}
//-----------------------------------------------------------------------------
void SparsityPattern::disp() const
{
  section("SparsityPattern");
  message("distributed : %u", distributed_);
  message("rank        : %u", rank_);
  begin(  "dim         :");
  for (uint i = 0; i < rank_; ++i)
  {
    message("%u : %u", i, dim_[i]);
  }
  end();
  message("range       : [ %8u, %8u [", local_range_[0], local_range_[1]);
  end();
}
//-----------------------------------------------------------------------------
void SparsityPattern::apply()
{
  finalized_ = true;

  if (!distributed_) return;

#ifdef HAVE_MPI

  message(1, "SparsityPattern : apply with %u remote rows", r_entries_.size());
  tic();

  /// Collect entries per owner
  uint const rank = MPI::rank();
  uint const pe_size = MPI::size();
  Array<uint> * sendbuf = new Array<uint>[pe_size];
  uint owner = 0;
  uint sendmax = 0;
  for (_ordered_map<uint, _ordered_set<uint> >::const_iterator it = r_entries_.begin();
       it != r_entries_.end(); ++it)
  {
    // Increment owner when jumping to another range
    while (range_[0][owner + 1] <= it->first)
    {
      sendmax = std::max(sendmax, (uint) sendbuf[owner].size());
      ++owner;
    }
    dolfin_assert(range_[0][owner] <= it->first);
    dolfin_assert(it->first < range_[0][owner + 1]);

    // Data packet [ global index, number of entries, [ indices ] ]
    sendbuf[owner].push_back(it->first);
    sendbuf[owner].push_back(it->second.size());
    for (_ordered_set<uint>::const_iterator c = it->second.begin();
         c != it->second.end(); ++c)
    {
      sendbuf[owner].push_back(*c);
    }
  }
  sendmax = std::max(sendmax, (uint) sendbuf[owner].size());
  uint recvmax = 0;
  MPI::all_reduce<MPI::max>(sendmax, recvmax);
  if (recvmax == 0)
  {
    return;
  }

  ///
  uint const r0 = local_range_[0][0];
  uint const r1 = local_range_[0][1];
  uint const c0 = local_range_[1][0];
  uint const c1 = local_range_[1][1];

  /// Exchange entries and add to diagonal and off-diagonal data structures
  uint * recvbuf = new uint[recvmax];
  for (uint j = 1; j < pe_size; ++j)
  {
    uint src = (rank - j + pe_size) % pe_size;
    uint dst = (rank + j) % pe_size;

    int recv_count = MPI::sendrecv( &sendbuf[dst][0], sendbuf[dst].size(), dst,
                                    recvbuf, recvmax, src, 1 );

    for (int k = 0; k < recv_count;)
    {
      uint const irow = recvbuf[k++];
      if (irow < r0 || r1 <= irow)
      {
        error("SparsityPattern : invalid global row number received %u", irow);
      }
      uint const ncol = recvbuf[k++];
      for (uint p = 0; p < ncol; ++p, ++k)
      {
        if (c0 <= recvbuf[k] && recvbuf[k] < c1)
        {
          d_entries_[irow - r0].insert(recvbuf[k]);
          ++d_count_;
        }
        else
        {
          o_entries_[irow - r0].insert(recvbuf[k]);
          ++o_count_;
        }
      }
    }
  }
  delete[] recvbuf;
  delete[] sendbuf;

  tocd(1);

#endif
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
