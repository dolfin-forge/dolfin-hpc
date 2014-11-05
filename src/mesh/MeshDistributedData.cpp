// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurélien Larcher, 2014.
//
// First added:  2008-07-03
// Last changed: 2014-03-18

#include <dolfin/config/dolfin_config.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Face.h>
#include <dolfin/main/MPI.h>
#include <dolfin/log/log.h>
#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshDistributedData::MeshDistributedData(MeshTopology& topology) :
    topology_(topology),
    topological_dim_(0),
    cell_dim_(0),
    facet_dim_(0),
    max_global_index_(0)
{
  for(uint i = 0; i < MAX_SIZE; ++i)
  {
    cached_global_indices_[i] = NULL;
  }
  set_invalid_numbering();
  set_invalid_ownership();
  clear();
  init(topology.dim());  // Set to zero is the mesh topology is uninitialized
}
//-----------------------------------------------------------------------------
MeshDistributedData::~MeshDistributedData()
{
  clear();
}
//-----------------------------------------------------------------------------
MeshDistributedData const& MeshDistributedData::operator=(
    MeshDistributedData const& other)
{
  clear();

  // Topological dimensions
  topological_dim_ = other.topological_dim_;
  cell_dim_ = other.cell_dim_;
  facet_dim_ = other.facet_dim_;

  //
  max_global_index_ = other.max_global_index_;
  for (uint i = 0; i <= cell_dim_; ++i)
  {
    // Numbering
    valid_numbering_[i] = other.valid_numbering_[i];
    num_global_[i] = other.num_global_[i];
    global_indices_[i] = other.global_indices_[i];
    local_indices_[i] = other.local_indices_[i];

    // Cached arrays
    finalized_[i] = other.finalized_[i];
    cached_global_size_[i] = other.cached_global_size_[i];
    if (finalized_[i])
    {
      dolfin_assert(cached_global_size_[i] > 0);
      cached_global_indices_[i] = new uint[cached_global_size_[i]];
      std::memcpy(cached_global_indices_[i], other.cached_global_indices_[i],
                  cached_global_size_[i] * sizeof(uint));
    }
    else
    {
      dolfin_assert(cached_global_size_[i] == 0);
      dolfin_assert(other.cached_global_indices_[i] == NULL);
      cached_global_indices_[i] = NULL;
    }

    // Ownership
    valid_ownership_[i] = other.valid_ownership_[i];
    adjacent_ranks_[i] = other.adjacent_ranks_[i];
    shared_[i] = other.shared_[i];
    shared_adj_[i] = other.shared_adj_[i];
    ghost_[i] = other.ghost_[i];
    ghost_owner_[i] = other.ghost_owner_[i];

    // Mapping
    valid_mapping_[i] = other.valid_mapping_[i];
    shared_mapping_[i] = other.shared_mapping_[i];
    ghost_mapping_[i] = other.ghost_mapping_[i];
  }

  return *this;
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::empty() const
{
  return (
      finalized_[0] ?
          (cached_global_size_[0] == 0) : (global_indices_[0].size() == 0));
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::is_finalized(uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return finalized_[dim];
}
//-----------------------------------------------------------------------------
void MeshDistributedData::init(uint dim)
{
  // Do not initialize at construction and if the mesh is a point cell
  if (dim > 0 && topological_dim_ == 0)
  {
    topological_dim_ = dim;
    cell_dim_ = topological_dim_;
    facet_dim_ = topological_dim_ - 1;
  }
  else if (dim != topological_dim_)
  {
    error("Setting a different topological dimension in distributed data.");
  }
}
//-----------------------------------------------------------------------------
void MeshDistributedData::clear()
{
  // Clear all numbering and ownership data structures
  for (uint i = 0; i < MAX_SIZE; ++i)
  {
    // Flushes numbering data structures and calls flush_numbering_cache
    flush_numbering_data(i);

    // Flushes ownership data structures and calls flush_mapping
    flush_ownership_data(i);
  }
  max_global_index_ = 0;
  topological_dim_ = 0;
  cell_dim_ = 0;
  facet_dim_ = 0;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::flush_numbering_data(uint dim)
{
  dolfin_assert(dim < MAX_SIZE);
  global_indices_[dim].clear();
  local_indices_[dim].clear();
  valid_numbering_[dim] = false;

  flush_numbering_cache(dim);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::flush_numbering_cache(uint dim)
{
  dolfin_assert(dim < MAX_SIZE);
  cached_global_size_[dim] = 0;
  if(cached_global_indices_[dim] != NULL)
  {
    dolfin_assert(finalized_[dim]);
    delete[] cached_global_indices_[dim];
  }
  cached_global_indices_[dim] = NULL;
  finalized_[dim] = false;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::flush_ownership_data(uint dim)
{
  dolfin_assert(dim < MAX_SIZE);
  adjacent_ranks_[dim].clear();
  shared_[dim].clear();
  shared_adj_[dim].clear();
  ghost_[dim].clear();
  ghost_owner_[dim].clear();
  valid_ownership_[dim] = false;

  flush_mapping(dim);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::flush_mapping(uint dim)
{
  dolfin_assert(dim < MAX_SIZE);
  shared_mapping_[dim].clear();
  ghost_mapping_[dim].clear();
  valid_mapping_[dim] = false;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::finalize(uint dim)
{
  if (dim > cell_dim_ || global_indices_[dim].size() == 0)
  {
    return;
  }
  if(finalized_[dim] == true)
  {
    error("Finalizing already finalized numbering for dimension %d", dim);
  }
  cached_global_size_[dim] = global_indices_[dim].size();
  dolfin_assert(cached_global_indices_[dim] == NULL);
  cached_global_indices_[dim] = new uint[global_indices_[dim].size()];
  for (_map<uint, uint>::iterator it = global_indices_[dim].begin();
  it != global_indices_[dim].end(); ++it)
  {
    cached_global_indices_[dim][it->first] = it->second;
  }
  max_global_index_ = std::max(max_global_index_, cached_global_size_[dim]);
  global_indices_[dim].clear();
  finalized_[dim] = true;
}
//-----------------------------------------------------------------------------
// Distributed entities numbering
//-----------------------------------------------------------------------------
bool MeshDistributedData::has_valid_numbering(uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return valid_numbering_[dim];
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_invalid_numbering()
{
  for (uint i = 0; i <= cell_dim_; ++i)
  {
    valid_numbering_[i] = false;
    finalized_[i] = false;
  }
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_map(uint local_index, uint global_index, uint dim)
{
  dolfin_assert(dim <= cell_dim_);
  global_indices_[dim][local_index] = global_index;
  local_indices_[dim][global_index] = local_index;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::apply_numbering(uint dim, _map<uint, uint> const& local,
                                        _map<uint, uint> const& global)
{
  dolfin_assert(dim <= cell_dim_);
  local_indices_[dim] = local;
  global_indices_[dim] = global;
  max_global_index_ = std::max(max_global_index_, num_global_[dim]);
  valid_numbering_[dim] = true;
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::num_global(uint dim) const
{
  if (dim > cell_dim_)
  {
    error("Trying to get global number of entities for invalid dimension.");
  }
  return num_global_[dim];
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::has_global(uint i, uint dim) const
{
  return (
      MPI::numProcesses() > 1 ?
          (local_indices_[dim].count(i) > 0) : (i < topology_.size(dim)));
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::has_global(MeshEntity const& entity) const
{
  return has_global(entity.index(), entity.dim());
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_global(uint i, uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  if (MPI::numProcesses() == 1)
  {
    return i;
  }
  else if (finalized_[dim])
  {
    dolfin_assert(cached_global_indices_[dim] != NULL);
    dolfin_assert(i < cached_global_size_[dim]);
    return cached_global_indices_[dim][i];
  }
  else
  {
    dolfin_assert(global_indices_[dim].count(i));
    return global_indices_[dim][i];
  }
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_global(MeshEntity const& e) const
{
  return get_global(e.index(), e.dim());
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_vertex_global(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }
  else if (finalized_[0])
  {
    dolfin_assert(cached_global_indices_[0] != NULL);
    dolfin_assert(i < cached_global_size_[0]);
    return cached_global_indices_[0][i];
  }
  else
  {
    dolfin_assert(global_indices_[0].count(i));
    return global_indices_[0][i];
  }
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_facet_global(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }
  else if (finalized_[facet_dim_])
  {
    dolfin_assert(cached_global_indices_[facet_dim_] != NULL);
    dolfin_assert(i < cached_global_size_[facet_dim_]);
    return cached_global_indices_[facet_dim_][i];
  }
  else
  {
    dolfin_assert(global_indices_[facet_dim_].count(i));
    return global_indices_[facet_dim_][i];
  }
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_cell_global(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }
  else if (finalized_[cell_dim_])
  {
    dolfin_assert(cached_global_indices_[cell_dim_] != NULL);
    dolfin_assert(i < cached_global_size_[cell_dim_]);
    return cached_global_indices_[cell_dim_][i];
  }
  else
  {
    dolfin_assert(global_indices_[cell_dim_].count(i));
    return global_indices_[cell_dim_][i];
  }
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_num_global(uint dim, uint const num_global)
{
  if (dim > cell_dim_)
  {
    error("Trying to set global number of entities for invalid dimension.");
  }
  if (num_global < topology_.size(dim))
  {
    error("Trying to set global number of entities lower than local number.");
  }
  num_global_[dim] = num_global;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::apply_num_global(uint dim, uint& offset)
{
  if (dim > cell_dim_)
  {
    error("Trying to set global number of entities for invalid dimension.");
  }
  if (MPI::numProcesses() > 1)
  {
#if HAVE_MPI
    offset = 0;
    uint num_owned = this->num_owned(dim);

#if ( MPI_VERSION > 1 )
    MPI_Exscan(&num_owned, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
    MPI_Scan(&num_owned, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
    offset -= num_owned;
#endif

    uint num_glb;
    MPI_Allreduce(&num_owned, &num_glb, 1, MPI_UNSIGNED, MPI_SUM,
                  MPI::DOLFIN_COMM);
    num_global_[dim] = num_glb;
  }
#endif
  else
  {
    offset = 0;
    num_global_[dim] = topology_.size(dim);
  }
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::has_local(uint i, uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return (
      MPI::numProcesses() > 1 ?
          (global_indices_[dim].count(i) > 0) : (i < topology_.size(dim)));
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::has_local(MeshEntity const& entity) const
{
  return has_local(entity.index(), entity.dim());
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_local(uint i, uint dim) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }
  dolfin_assert(dim <= cell_dim_);
  dolfin_assert(local_indices_[dim].count(i));
  return local_indices_[dim][i];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_local(MeshEntity const& e) const
{
  return get_local(e.index(), e.dim());
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_vertex_local(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }
  dolfin_assert(local_indices_[0].count(i));
  return local_indices_[0][i];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_facet_local(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  dolfin_assert(local_indices_[facet_dim_].count(i));
  return local_indices_[facet_dim_][i];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_cell_local(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  dolfin_assert(local_indices_[cell_dim_].count(i));
  return local_indices_[cell_dim_][i];
}
//-----------------------------------------------------------------------------
// Distributed entities ownership
//-----------------------------------------------------------------------------
bool MeshDistributedData::has_valid_ownership(uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return valid_ownership_[dim];
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::has_valid_mapping(uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return valid_mapping_[dim];
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_invalid_ownership()
{
  for (uint i = 0; i < cell_dim_; ++i)
  {
    valid_ownership_[i] = false;
  }
}
//-----------------------------------------------------------------------------
void MeshDistributedData::apply_ownership(uint dim)
{
  dolfin_assert(dim <= cell_dim_);
  valid_ownership_[dim] = true;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::remap_ownership(int const* mapping)
{
  for (uint i = 0; i < cell_dim_; ++i)
  {
    for (MeshGhostIterator it(*this, i); !it.end(); ++it)
    {
      set_ghost_owner(it.index(), mapping[it.owner()], i);
    }
    for (MeshSharedIterator it(*this, i); !it.end(); ++it)
    {
      _set<uint> new_adj;
      _set<uint> const& adj = it.adj();
      for(_set<uint>::const_iterator adjit = adj.begin(); adjit != adj.end();
          ++adjit)
      {
        new_adj.insert(mapping[*adjit]);
      }
      setall_shared_adj(it.index(), new_adj, i);
    }
#ifdef ENABLE_P1_OPTIMIZATIONS
    break;
#endif
  }
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::is_shared(uint i, uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return (MPI::numProcesses() > 1 ? (shared_[dim].count(i) > 0) : false);
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::is_shared(MeshEntity const& entity) const
{
  return is_shared(entity.index(), entity.dim());
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::num_shared(uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return shared_[dim].size();
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::num_shared_with(uint rank, uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  AdjacentMapping::const_iterator it = shared_mapping_[dim].find(rank);
  if (it != shared_mapping_[dim].end())
  {
    return it->second.first.size();
  }
  return 0;
}
//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_shared_mapping_to(uint rank,
                                                              uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  dolfin_assert(adjacent_ranks_[dim].count(rank) > 0);
  dolfin_assert(shared_mapping_[dim].find(rank)->second.first.size()
                    == shared_mapping_[dim].find(rank)->second.second.size());
  AdjacentMapping::const_iterator it = shared_mapping_[dim].find(rank);
  if (it == shared_mapping_[dim].end())
  {
    if(adjacent_ranks_[dim].count(rank) > 0)
    {
      error("Shared mapping does not exists for adjacent rank %d", rank);
    }
    else
    {
      error("Requesting shared mapping from non-adjacent rank %d", rank);
    }
  }
  return it->second.first;
}
//-----------------------------------------------------------------------------
Array<uint>& MeshDistributedData::get_shared_mapping_to(uint rank, uint dim)
{
  dolfin_assert(dim <= cell_dim_);
  if(adjacent_ranks_[dim].count(rank) == 0)
  {
    error("Requesting shared mapping from non-adjacent rank %d", rank);
  }
  return shared_mapping_[dim][rank].first;
}
//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_shared_mapping_from(uint rank,
                                                                uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  dolfin_assert(adjacent_ranks_[dim].count(rank) > 0);
  dolfin_assert(shared_mapping_[dim].find(rank)->second.first.size()
                  == shared_mapping_[dim].find(rank)->second.second.size());
  AdjacentMapping::const_iterator it = shared_mapping_[dim].find(rank);
  if (it == shared_mapping_[dim].end())
  {
    if(adjacent_ranks_[dim].count(rank) > 0)
    {
      error("Shared mapping does not exists for adjacent rank %d", rank);
    }
    else
    {
      error("Requesting shared mapping from non-adjacent rank %d", rank);
    }
  }
  return it->second.second;
}
//-----------------------------------------------------------------------------
Array<uint>& MeshDistributedData::get_shared_mapping_from(uint rank, uint dim)
{
  dolfin_assert(dim <= cell_dim_);
  if(adjacent_ranks_[dim].count(rank) == 0)
  {
    error("Requesting shared mapping from non-adjacent rank %d", rank);
  }
  return shared_mapping_[dim][rank].second;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared(MeshEntity const& m)
{
  set_shared(m.index(), m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared(uint local_index, uint dim)
{
  dolfin_assert(dim <= cell_dim_);
  shared_[dim].insert(local_index);
}
//-----------------------------------------------------------------------------
_set<uint> const& MeshDistributedData::get_adj_ranks(uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return adjacent_ranks_[dim];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::num_adj_ranks(uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return adjacent_ranks_[dim].size();
}
//-----------------------------------------------------------------------------
_set<uint> const& MeshDistributedData::get_shared_adj(MeshEntity const& m) const
{
  return get_shared_adj(m.index(), m.dim());
}
//-----------------------------------------------------------------------------
_set<uint> const& MeshDistributedData::get_shared_adj(uint local_index,
    uint dim) const
{
  dolfin_assert(is_shared(local_index, dim));
  return shared_adj_[dim][local_index];
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared_adj(MeshEntity const& m, uint rank)
{
  set_shared_adj(m.index(), rank, m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared_adj(uint i, uint rank, uint dim)
{
  dolfin_assert(dim <= cell_dim_);
  shared_adj_[dim][i].insert(rank);
  adjacent_ranks_[dim].insert(rank);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::setall_shared_adj(uint i, _set<uint> const& ranks,
uint dim)
{
  dolfin_assert(dim <= cell_dim_);
  shared_adj_[dim][i].clear();
  shared_adj_[dim][i].insert(ranks.begin(), ranks.end());
  adjacent_ranks_[dim].insert(ranks.begin(), ranks.end());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::setall_shared_adj(MeshEntity const& m,
                                            _set<uint> const& ranks)
{
  setall_shared_adj(m.index(), ranks, m.dim());
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::is_ghost(uint i, uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return (MPI::numProcesses() > 1 ? (ghost_[dim].count(i) > 0) : false);
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::is_ghost(MeshEntity const& entity) const
{
  return is_ghost(entity.index(), entity.dim());
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::num_owned(uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return topology_.size(dim) - ghost_[dim].size();
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::num_ghost(uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  return ghost_[dim].size();
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::num_ghost_from(uint rank, uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  AdjacentMapping::const_iterator it = ghost_mapping_[dim].find(rank);
  if (it != ghost_mapping_[dim].end())
  {
    return it->second.first.size();
  }
  return 0;
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_owner(uint local_index, uint dim) const
{
  if (MPI::numProcesses() == 1)
  {
    return 0;
  }
  dolfin_assert(dim <= cell_dim_);
  dolfin_assert(ghost_owner_[dim].count(local_index));
  return ghost_owner_[dim][local_index];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_owner(MeshEntity const& e) const
{
  return get_owner(e.index(), e.dim());
}
//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_ghost_mapping_to(uint rank,
                                                             uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  dolfin_assert(adjacent_ranks_[dim].count(rank) > 0);
  dolfin_assert(ghost_mapping_[dim].find(rank)->second.first.size()
                  == ghost_mapping_[dim].find(rank)->second.second.size());
  AdjacentMapping::const_iterator it = ghost_mapping_[dim].find(rank);
  if (it == ghost_mapping_[dim].end())
  {
    if(adjacent_ranks_[dim].count(rank) > 0)
    {
      error("Ghost mapping does not exists for adjacent rank %d", rank);
    }
    else
    {
      error("Requesting ghost mapping from non-adjacent rank %d", rank);
    }
  }
  return it->second.first;
}
//-----------------------------------------------------------------------------
Array<uint>& MeshDistributedData::get_ghost_mapping_to(uint rank, uint dim)
{
  dolfin_assert(dim <= cell_dim_);
  if(adjacent_ranks_[dim].count(rank) == 0)
  {
    error("Requesting ghost mapping from non-adjacent rank %d", rank);
  }
  return ghost_mapping_[dim][rank].first;
}
//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_ghost_mapping_from(uint rank,
                                                               uint dim) const
{
  dolfin_assert(dim <= cell_dim_);
  dolfin_assert(adjacent_ranks_[dim].count(rank) > 0);
  dolfin_assert(ghost_mapping_[dim].find(rank)->second.first.size()
                == ghost_mapping_[dim].find(rank)->second.second.size());
  AdjacentMapping::const_iterator it = ghost_mapping_[dim].find(rank);
  if (it == ghost_mapping_[dim].end())
  {
    if(adjacent_ranks_[dim].count(rank) > 0)
    {
      error("Ghost mapping does not exists for adjacent rank %d", rank);
    }
    else
    {
      error("Requesting ghost mapping from non-adjacent rank %d", rank);
    }
  }
  return it->second.second;
}
//-----------------------------------------------------------------------------
Array<uint>& MeshDistributedData::get_ghost_mapping_from(uint rank, uint dim)
{
  dolfin_assert(dim <= cell_dim_);
  if(adjacent_ranks_[dim].count(rank) == 0)
  {
    error("Requesting ghost mapping from non-adjacent rank %d", rank);
  }
  return ghost_mapping_[dim][rank].second;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost(MeshEntity const& m)
{
  set_ghost(m.index(), m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost(uint local_index, uint dim)
{
  dolfin_assert(dim <= cell_dim_);
  set_shared(local_index, dim);
  ghost_[dim].insert(local_index);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost_owner(MeshEntity const& m, uint rank)
{
  set_ghost_owner(m.index(), rank, m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost_owner(uint i, uint rank, uint dim)
{
  dolfin_assert(dim <= cell_dim_);
  set_shared_adj(i, rank, dim);
  ghost_owner_[dim][i] = rank;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::disp() const
{
  cout << "MeshDistributedData" << endl;
  cout << "-------------------" << endl;

  begin("");
  cout << "Topological dimension     : " << (uint) topological_dim_ << endl;
  cout << "Cell dimension            : " << (uint) cell_dim_ << endl;
  cout << "Facet dimension           : " << (uint) facet_dim_ << endl;
  skip();
  cout << "Maximum global index      : " << (uint) max_global_index_ << endl;
  cout << "Number of global vertices : " << (uint) num_global_[0] << endl;
  cout << "Valid vertex numbering    : " << (bool) valid_numbering_[0] << endl;
  if (topological_dim_ > 1)
  {
    cout << "Number of global edges    : " << (uint) num_global_[1] << endl;
    cout << "Valid edge   numbering    : " << (bool) valid_numbering_[1]
         << endl;
  }
  if (topological_dim_ > 2)
  {
    cout << "Number of global faces    : " << (uint) num_global_[2] << endl;
    cout << "Valid face   numbering    : " << (bool) valid_numbering_[2]
         << endl;
  }
  cout << "Number of global cells    : " << (uint) num_global_[cell_dim_]
       << endl;
  cout << "Valid cell   numbering    : " << (bool) valid_numbering_[cell_dim_]
       << endl;
  skip();
  cout << "Number of shared entities : " << endl;
  for (uint d = 0; d < topological_dim_; ++d)
  {
    cout << "  - dim " << " : " << (uint) this->num_shared(d) << endl;
  }
  skip();
  cout << "Number of ghost entities : " << endl;
  for (uint d = 0; d < topological_dim_; ++d)
  {
    cout << "  - dim " << " : " << (uint) this->num_ghost(d) << endl;
  }
  end();
}

}

