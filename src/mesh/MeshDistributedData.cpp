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
#include <string.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshDistributedData::MeshDistributedData(MeshTopology& topology) :
    topological_dim_(0),
    cell_dim_(0),
    facet_dim_(0),
    max_global_vertex_index_(0),
    valid_edge_ownership_(false),
    valid_face_ownership_(false),
    valid_shared_facets_mapping_(false),
    finalized_(false),
    global_vertex_indices_(0),
    global_facet_indices_(0),
    global_cell_indices_(0),
    global_vertex_indices_size_(0),
    global_facet_indices_size_(0),
    global_cell_indices_size_(0)
{
  init(topology.dim()); // Set to zero is the mesh topology is uninitialized
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

  topological_dim_ = other.topological_dim_;
  cell_dim_ = other.cell_dim_;
  facet_dim_ = other.facet_dim_;

  max_global_vertex_index_ = other.max_global_vertex_index_;

  valid_edge_ownership_ = other.valid_edge_ownership_;
  valid_face_ownership_ = other.valid_face_ownership_;

  valid_shared_facets_mapping_ = other.valid_shared_facets_mapping_;

  for (uint i = 0; i < MAX_DIM + 1; ++i)
  {
    num_global_[i] = other.num_global_[i];
    valid_numbering_[i] = other.valid_numbering_[i];
    global_indices_[i] = other.global_indices_[i];
    local_indices_[i] = other.local_indices_[i];
  }

  for (uint i = 0; i < MAX_DIM; ++i)
  {
    adjacent_ranks_[i] = other.adjacent_ranks_[i];
    shared_[i] = other.shared_[i];
    shared_adj_[i] = other.shared_adj_[i];
    ghost_[i] = other.ghost_[i];
    ghost_owner_[i] = other.ghost_owner_[i];
    shared_mapping_[i] = other.shared_mapping_[i];
    ghost_mapping_[i] = other.ghost_mapping_[i];
  }

  finalized_ = other.finalized_;

  global_vertex_indices_size_ = other.global_vertex_indices_size_;
  global_facet_indices_size_ = other.global_facet_indices_size_;
  global_cell_indices_size_ = other.global_cell_indices_size_;

  if (finalized_)
  {
    dolfin_assert(global_vertex_indices_size_ > 0);
    dolfin_assert(global_cell_indices_size_ > 0);

    global_vertex_indices_ = new uint[global_vertex_indices_size_];
    memcpy(global_vertex_indices_, other.global_vertex_indices_,
           global_vertex_indices_size_ * sizeof(uint));

    global_facet_indices_ = new uint[global_facet_indices_size_];
    memcpy(global_facet_indices_, other.global_facet_indices_,
           global_facet_indices_size_ * sizeof(uint));

    global_cell_indices_ = new uint[global_cell_indices_size_];
    memcpy(global_cell_indices_, other.global_cell_indices_,
           global_cell_indices_size_ * sizeof(uint));
  }

  return *this;
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::empty() const
{
  return (global_indices_[0].size() == 0);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::init(uint const dim)
{
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
  topological_dim_ = 0;
  cell_dim_ = 0;
  facet_dim_ = 0;
  max_global_vertex_index_ = 0;

  for (uint i = 0; i < MAX_DIM; ++i)
  {
    adjacent_ranks_[i].clear();
    shared_[i].clear();
    shared_adj_[i].clear();
    ghost_[i].clear();
    ghost_owner_[i].clear();
    shared_mapping_[i].clear();
    ghost_mapping_[i].clear();
  }

  for (uint i = 0; i < MAX_DIM + 1; ++i)
  {
    global_indices_[i].clear();
    local_indices_[i].clear();
    valid_numbering_[i] = false;
  }

  valid_edge_ownership_ = false;
  valid_face_ownership_ = false;

  valid_shared_facets_mapping_ = false;

  delete[] global_vertex_indices_;
  global_vertex_indices_ = NULL;

  delete[] global_facet_indices_;
  global_facet_indices_ = NULL;

  delete[] global_cell_indices_;
  global_cell_indices_ = NULL;

  finalized_ = false;

}
//-----------------------------------------------------------------------------
void MeshDistributedData::finalize(uint const dim)
{

  _map<uint, uint>::iterator it;

  if (dim == 0) // Vertices
  {
    delete[] global_vertex_indices_;
    global_vertex_indices_ = new uint[global_indices_[0].size()];

    for(it = global_indices_[0].begin(); it != global_indices_[0].end(); ++it)
    {
      global_vertex_indices_[it->first] = it->second;
    }
    global_vertex_indices_size_ = global_indices_[0].size();
    max_global_vertex_index_ = global_vertex_indices_size_;
    global_indices_[0].clear();
  }
  else if (dim == facet_dim_) // Facets
  {
    delete[] global_facet_indices_;
    global_facet_indices_ = new uint[global_indices_[dim-1].size()];

    for(it = global_indices_[dim-1].begin(); it != global_indices_[dim-1].end(); ++it)
    {
      global_facet_indices_[it->first] = it->second;
    }
    global_facet_indices_size_ = global_indices_[2].size();
    max_global_vertex_index_ = global_facet_indices_size_;
    global_indices_[2].clear();
  }
  else if (dim == cell_dim_) // Cells
  {
    delete[] global_cell_indices_;
    global_cell_indices_ = new uint[global_indices_[dim].size()];

    for(it = global_indices_[dim].begin(); it != global_indices_[dim].end(); ++it)
    {
      global_cell_indices_[it->first] = it->second;
    }
    global_cell_indices_size_ = global_indices_[MAX_DIM].size();
    global_indices_[MAX_DIM].clear();
  }
  else
  {
    error("MeshDistributedData::finalize not implemented for %ud.", dim);
  }

  finalized_ = true;
}
//-----------------------------------------------------------------------------
// Distributed entities numbering
//-----------------------------------------------------------------------------
bool MeshDistributedData::has_global(uint i, uint dim) const
{
  return (MPI::numProcesses() > 1 ? (local_indices_[dim].count(i) > 0) : true);
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::has_global(MeshEntity const& entity) const
{
  return has_global(entity.index(), entity.dim());
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_global(uint i, uint dim) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  if (dim == 0 && finalized_)
  {
    return global_vertex_indices_[i];
  }
  else
  {
    dolfin_assert( global_indices_[dim].count(i) );
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

  if (finalized_)
  {
    return global_vertex_indices_[i];
  }
  else
  {
    dolfin_assert( global_indices_[0].count(i) );
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

  if (finalized_)
  {
    return global_facet_indices_[i];
  }
  else
  {
    dolfin_assert( facet_dim_ != 0 );dolfin_assert( global_indices_[facet_dim_].count(i) );
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

  if (finalized_)
  {
    return global_cell_indices_[i];
  }
  else
  {
    dolfin_assert( cell_dim_ != 0 );dolfin_assert( global_indices_[cell_dim_].count(i) );
    return global_indices_[cell_dim_][i];
  }

}
//-----------------------------------------------------------------------------
bool MeshDistributedData::has_local(uint i, uint dim) const
{
  return (MPI::numProcesses() > 1 ? (global_indices_[dim].count(i) > 0) : true);
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

  dolfin_assert( local_indices_[dim].count(i) );
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

  dolfin_assert( local_indices_[0].count(i) );
  return local_indices_[0][i];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_facet_local(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  dolfin_assert( facet_dim_ != 0 );dolfin_assert( local_indices_[facet_dim_].count(i) );
  return local_indices_[facet_dim_][i];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_cell_local(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  dolfin_assert( cell_dim_ != 0 );dolfin_assert( local_indices_[cell_dim_].count(i) );
  return local_indices_[cell_dim_][i];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::num_global(uint dim) const
{
  if (dim > topological_dim_)
  {
    error("Trying to set global number of entities for invalid dimension.");
  }
  return num_global_[dim];
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_map(uint local_index, uint global_index, uint dim)
{
  if (dim == 0)
  {
    max_global_vertex_index_ = std::max(max_global_vertex_index_, global_index);
  }
  global_indices_[dim][local_index] = global_index;
  local_indices_[dim][global_index] = local_index;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_num_global(uint dim, uint num_global)
{
  if (dim > topological_dim_)
  {
    error("Trying to set global number of entities for invalid dimension.");
  }
  num_global_[dim] = num_global;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_invalid_numbering()
{
  std::fill(&valid_numbering_[0], &valid_numbering_[MAX_DIM], false);
  finalized_ = false;
}
//-----------------------------------------------------------------------------
// Distributed entities ownership
//-----------------------------------------------------------------------------
_set<uint> const& MeshDistributedData::get_adj(uint dim) const
{
  return adjacent_ranks_[dim];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::num_adj(uint dim) const
{
  return adjacent_ranks_[dim].size();
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::is_shared(uint i, uint dim) const
{
  return (MPI::numProcesses() > 1 ? (shared_[dim].count(i) > 0) : false);
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::num_shared(uint dim) const
{
  return shared_[dim].size();
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared(MeshEntity const& m)
{
  set_shared(m.index(), m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared(uint local_index, uint dim)
{
  shared_[dim].insert(local_index);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost(MeshEntity const& m)
{
  set_ghost(m.index(), m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost(uint local_index, uint dim)
{
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
  set_shared_adj(i, rank, dim);
  ghost_owner_[dim][i] = rank;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared_adj(MeshEntity const& m, uint rank)
{
  set_shared_adj(m.index(), rank, m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared_adj(uint i, uint rank, uint dim)
{
  shared_adj_[dim][i].insert(rank);
  adjacent_ranks_[dim].insert(rank);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::setall_shared_adj(uint i, _set<uint> const& ranks,
uint dim)
{
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
uint MeshDistributedData::get_owner(MeshEntity const& e) const
{
  return get_owner(e.index(), e.dim());
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_owner(uint local_index, uint dim) const
{
  if (MPI::numProcesses() == 1)
  {
    return 0;
  }dolfin_assert( ghost_owner_[dim].count(local_index) );
  return ghost_owner_[dim][local_index];
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
uint MeshDistributedData::num_shared_with(uint rank, uint dim) const
{
  AdjacentMapping::const_iterator it = shared_mapping_[dim].find(rank);
  if (it != shared_mapping_[dim].end())
  {
    return it->second.first.size();
  }
  return 0;
}

//-----------------------------------------------------------------------------
uint MeshDistributedData::num_ghost_from(uint rank, uint dim) const
{
  AdjacentMapping::const_iterator it = ghost_mapping_[dim].find(rank);
  if (it != ghost_mapping_[dim].end())
  {
    return it->second.first.size();
  }
  return 0;
}

//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_shared_mapping_to(uint rank,
                                                              uint dim) const
{
  dolfin_assert(shared_mapping_[dim].find(rank) != shared_mapping_[dim].end());
  return shared_mapping_[dim].find(rank)->second.first;
}

//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_shared_mapping_from(uint rank,
                                                                uint dim) const
{
  dolfin_assert(shared_mapping_[dim].find(rank) != shared_mapping_[dim].end());
  return shared_mapping_[dim].find(rank)->second.second;
}

//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_ghost_mapping_to(uint rank,
                                                             uint dim) const
{
  dolfin_assert(shared_mapping_[dim].find(rank) != shared_mapping_[dim].end());
  return ghost_mapping_[dim].find(rank)->second.first;
}

//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_ghost_mapping_from(uint rank,
                                                               uint dim) const
{
  dolfin_assert(shared_mapping_[dim].find(rank) != shared_mapping_[dim].end());
  return ghost_mapping_[dim].find(rank)->second.second;
}

//-----------------------------------------------------------------------------
void MeshDistributedData::remap_owner(int* mapping)
{

  for (uint i = 0; i < MAX_DIM; i++)
  {
    for (MeshGhostIterator it(*this, i); !it.end(); ++it)
    {
      //FIXME: logic with shared_adj
      set_ghost_owner(it.index(), mapping[it.owner()], i);
    }
#ifdef ENABLE_P1_OPTIMIZATIONS
    break;
#endif
  }

}
//-----------------------------------------------------------------------------
bool MeshDistributedData::is_shared(MeshEntity const& entity) const
{
  return is_shared(entity.index(), entity.dim());
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::is_ghost(MeshEntity const& entity) const
{
  return is_ghost(entity.index(), entity.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::flush_mappings(uint dim)
{
  shared_mapping_[dim].clear();
  ghost_mapping_[dim].clear();
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
  cout << "Maximum global index      : " << (uint) max_global_vertex_index_
       << endl;
  cout << "Number of global vertices : " << (uint) num_global_[0] << endl;
  cout << "Number of global edges    : " << (uint) num_global_[1] << endl;
  cout << "Number of global faces    : " << (uint) num_global_[2] << endl;
  cout << "Number of global cells    : " << (uint) num_global_[cell_dim_]
       << endl;
  skip();
  cout << "Valid vertex numbering    : " << (bool) valid_numbering_[0] << endl;
  cout << "Valid edge   numbering    : " << (bool) valid_numbering_[1] << endl;
  cout << "Valid face   numbering    : " << (bool) valid_numbering_[2] << endl;
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

