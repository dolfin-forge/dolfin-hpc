// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Magnus Vikstrom, 2008.
// Modified by Anders Logg, 2008.
// Modified by Niclas Jansson, 2009.
// Modified by Aurélien Larcher, 2015.
//
// First added:  2007-03-13
// Last changed: 2009-03-03

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/main/MPI.h>
#include <iostream>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
SparsityPattern::SparsityPattern(uint rank, uint const * dims,
                                 bool distributed) :
    range(NULL),
    initialized_(false),
    finalized_(false),
    blocked_(false),
    distributed_(distributed)
{
  init(rank, dims, distributed);
}
//-----------------------------------------------------------------------------
SparsityPattern::SparsityPattern() :
    range(NULL),
    initialized_(false),
    finalized_(false),
    blocked_(false),
    distributed_(false)
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
  initialized_ = false;
  finalized_ = false;
  blocked_ = false;
  distributed_ = false;
  delete[] range;
  range = NULL;
  dim[0] = 0;
  dim[1] = 0;
  off_processor.clear();
  sparsity_pattern.clear();
  o_sparsity_pattern.clear();
}
//-----------------------------------------------------------------------------
void SparsityPattern::init(uint rank, uint const * dims, bool distributed)
{
  dolfin_assert(rank <= 2);
  if(initialized_)
  {
    error("Sparsity pattern has already been initialized");
  }
  distributed_ = distributed;
  for (uint i = 0; i < rank; ++i)
  {
    dim[i] = dims[i];
  }
  if(rank == 2 && distributed_ && dim[0] != dim[1])
  {
    error("Only square matrices are supported in parallel");
  }
  sparsity_pattern.clear();
  if(distributed_)
  {
    o_sparsity_pattern.clear();
  }

  // Initialize range
  uint num_procs = dolfin::MPI::numProcesses();
  range = new uint[num_procs + 1];
  range[0] = 0;
  for (uint p = 0; p < num_procs; ++p)
  {
    range[p + 1] = range[p] + dim[0] / num_procs
                    + ((dim[0] % num_procs) > p ? 1 : 0);
  }
  initialized_ = true;
}
//-----------------------------------------------------------------------------
void SparsityPattern::init(uint rank, uint const * dims)
{
  init(rank, dims, false);
}
//-----------------------------------------------------------------------------
void SparsityPattern::pinit(uint rank, uint const * dims)
{
  init(rank, dims, true);
}
//-----------------------------------------------------------------------------
void SparsityPattern::insert(uint const * num_rows, uint const * const * rows)
{
  dolfin_assert(!distributed_);
  for (uint i = 0; i < num_rows[0]; ++i)
  {
    for (uint j = 0; j < num_rows[1]; ++j)
    {
      sparsity_pattern[rows[0][i]].insert(rows[1][j]);
    }
  }
}
//-----------------------------------------------------------------------------
void SparsityPattern::pinsert(uint const * num_rows, uint const * const * rows)
{
  dolfin_assert(distributed_);
  uint process = dolfin::MPI::processNumber();

  for (uint i = 0; i < num_rows[0]; ++i)
  {
    uint const global_row = rows[0][i];
    // If not in a row "owned" by this processor
    if (global_row < range[process] || global_row >= range[process + 1])
    {

      off_processor.push_back(global_row);
      off_processor.push_back(num_rows[1]);
      for (uint j = 0; j < num_rows[1]; ++j)
      {
        off_processor.push_back(rows[1][j]);
      }

      continue;
    }

    for (uint j = 0; j < num_rows[1]; ++j)
    {
      uint const global_col = rows[1][j];
      // On the off-diagonal
      if (global_col < range[process] || global_col >= range[process + 1])
      {
        o_sparsity_pattern[rows[0][i]].insert(rows[1][j]);
      }
      // On the diagonal
      else
      {
        sparsity_pattern[rows[0][i]].insert(rows[1][j]);
      }
    }
  }

}
//-----------------------------------------------------------------------------
dolfin::uint SparsityPattern::size(uint n) const
{
  dolfin_assert(n < 2);
  if (distributed_)
  {
    return (range[MPI::processNumber() + 1] - range[MPI::processNumber()]);
  }
  else
  {
    return dim[n];
  }
}
//-----------------------------------------------------------------------------
void SparsityPattern::numNonZeroPerRow(uint nzrow[]) const
{
  if (dim[1] == 0)
  {
    error("Non-zero entries per row can be computed for matrices only.");
  }

  if (sparsity_pattern.size() == 0)
  {
    error("Sparsity pattern has not been computed.");
  }

  // Compute number of nonzeros per row
#if __SUNPRO_CC
  std::map<uint, std::set<int> >::const_iterator set;
#else
  std::map<uint const, std::set<int> >::const_iterator set;
#endif
  for (set = sparsity_pattern.begin(); set != sparsity_pattern.end(); ++set)
  {
    nzrow[set->first - sparsity_pattern.begin()->first] = set->second.size();
  }

}
//-----------------------------------------------------------------------------
void SparsityPattern::numNonZeroPerRow(uint process_number, uint d_nzrow[],
                                       uint o_nzrow[]) const
{
  if (dim[1] == 0)
  {
    error("Non-zero entries per row can be computed for matrices only.");
  }

  if (sparsity_pattern.size() == 0) error(
      "Sparsity pattern has not been computed.");
#if __SUNPRO_CC
  std::map<uint, std::set<int> >::const_iterator it;
#else
  std::map<uint const, std::set<int> >::const_iterator it;
#endif
  // Compute number of nonzeros per row diagonal and off-diagonal
  uint offset = range[process_number];
  for (uint i = 0; i + offset < range[process_number + 1]; ++i)
  {

    it = sparsity_pattern.find(i + offset);
    if (it != sparsity_pattern.end())
    {
      d_nzrow[i] = it->second.size();
    }
    else
    {
      d_nzrow[i] = 0;
    }
    it = o_sparsity_pattern.find(i + offset);
    if (it != o_sparsity_pattern.end())
    {
      o_nzrow[i] = it->second.size();
    }
    else
    {
      o_nzrow[i] = 0;
    }

    //    d_nzrow[i] = sparsity_pattern[i+offset].size();
    //    o_nzrow[i] = o_sparsity_pattern[i+offset].size();
  }
}
//-----------------------------------------------------------------------------
dolfin::uint SparsityPattern::numNonZero() const
{
  if (dim[1] == 0) error(
      "Total non-zeros entries can be computed for matrices only.");

  if (sparsity_pattern.size() == 0) error(
      "Sparsity pattern has not been computed.");

  // Compute total number of nonzeros per row
  uint nz = 0;
#if __SUNPRO_CC
  std::map<uint, std::set<int> >::const_iterator set;
#else
  std::map<uint const, std::set<int> >::const_iterator set;
#endif
  for (set = sparsity_pattern.begin(); set != sparsity_pattern.end(); ++set)
    nz += set->second.size();
  return nz;

  return 0;
}
//-----------------------------------------------------------------------------
void SparsityPattern::disp() const
{
  cout << "SparsityPattern" << endl;
  cout << "---------------" << endl;
  cout << endl;
  if (MPI::numProcesses() == 1)
  {
   if ( dim[1] == 0 )
   {
     warning("Only matrix sparsity patterns can be displayed.");
   }

#if __SUNPRO_CC
    std::map<uint,  std::set<int> >::const_iterator set;
#else
    std::map<uint const,  std::set<int> >::const_iterator set;
#endif
   std::set<int>::const_iterator jj;
   for(set = sparsity_pattern.begin(); set != sparsity_pattern.end(); ++set)
   {
     cout << set->first << "\t: (" << (uint) set->second.size() << ")";
     for(jj = set->second.begin(); jj != set->second.end(); ++jj)
     {
       cout << "\t" << *jj;
     }
     cout << endl;
   }
  }
}
//-----------------------------------------------------------------------------
void SparsityPattern::processRange(uint process_number, uint local_range[])
{
  local_range[0] = range[process_number];
  local_range[1] = range[process_number + 1];
}
//-----------------------------------------------------------------------------
dolfin::uint SparsityPattern::numLocalRows(uint process_number) const
{
  return range[process_number + 1] - range[process_number];
}
//-----------------------------------------------------------------------------
void SparsityPattern::initRange(uint num_local)
{
  if (range) delete[] range;

  if(distributed_)
  {
    uint num_procs = dolfin::MPI::numProcesses();
    range = new uint[num_procs + 1];
    range[0] = 0;

    uint *local = new uint[num_procs];
    local[dolfin::MPI::processNumber()] = num_local;

#ifdef HAVE_MPI
    MPI_Allgather(&num_local, 1, MPI_UNSIGNED, local, 1, MPI_UNSIGNED,
                  MPI::DOLFIN_COMM);
#endif

    for (uint p = 0; p < num_procs; ++p)
    {
      range[p + 1] = range[p] + local[p];
    }

    delete[] local;
  }
  else
  {
    range = new uint[2];
    range[0] = 0;
    range[1] = num_local;
  }
}
//-----------------------------------------------------------------------------
void SparsityPattern::apply()
{
  finalized_ = true;

  if (!distributed_) return;

#ifdef HAVE_MPI
  int rank = MPI::processNumber();
  int pe_size = MPI::numProcesses();

  MPI_Status status;
  int maxoff, recv_count, global_row, src, dest;

  int numoff = off_processor.size();
  MPI_Allreduce(&numoff, &maxoff, 1, MPI_INT, MPI_MAX, MPI::DOLFIN_COMM);
  int *recv_buff = new int[maxoff];

  for (int j = 1; j < pe_size; j++)
  {

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&off_processor[0], numoff, MPI_INT, dest, 1, recv_buff, maxoff,
                 MPI_INT, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int i = 0; i < recv_count;)
    {

      global_row = recv_buff[i];

      if (global_row >= (int) range[rank] && global_row < (int) range[rank + 1])
      {
        for (int k = 0; k < recv_buff[i + 1]; k++)
        {
          if (recv_buff[k + i + 2] < (int) range[rank]
              || recv_buff[k + i + 2] >= (int) range[rank + 1])
          {
            o_sparsity_pattern[global_row].insert(recv_buff[k + i + 2]);
          }
          else
          {
            sparsity_pattern[global_row].insert(recv_buff[k + i + 2]);
          }
        }

      }
      i += 2 + recv_buff[i + 1];
    }

  }

  delete[] recv_buff;
#endif
}
//-----------------------------------------------------------------------------
Array< _set<int> > const& SparsityPattern::pattern() const
{
  if (finalized_)
  {
    error("Sparsity pattern has not been finalized");
  }
  return pattern_;
}
//-----------------------------------------------------------------------------

}

