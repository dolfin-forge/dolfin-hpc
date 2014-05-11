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
    _dim(0),
    _cell_dim(0),
    _facet_dim(0),
    _max_global_index(0),
    _num_global_vertex(0),
    _num_global_edge(0),
    _num_global_face(0),
    _num_global_cell(0),
    _valid_vertex_numbering(false),
    _valid_edge_numbering(false),
    _valid_face_numbering(false),
    _valid_cell_numbering(false),
    _valid_edge_ownership(false),
    _valid_face_ownership(false),
    _valid_shared_facets_mapping(false),
    _finalized(false),
    _global_vertex_indices(0),
    _global_facet_indices(0),
    _global_cell_indices(0),
    _global_vertex_indices_size(0),
    _global_facet_indices_size(0),
    _global_cell_indices_size(0)
{
  init(topology.dim()); // Set to zero is the mesh topology is uninitialized
}
//-----------------------------------------------------------------------------
MeshDistributedData::~MeshDistributedData()
{
  clear();
}
//-----------------------------------------------------------------------------
const MeshDistributedData& MeshDistributedData::operator=(
    const MeshDistributedData& distributed_data)
{
  clear();

  _dim = distributed_data._dim;
  _cell_dim = distributed_data._cell_dim;
  _facet_dim = distributed_data._facet_dim;

  _max_global_index = distributed_data._max_global_index;

  _valid_vertex_numbering = distributed_data._valid_vertex_numbering;
  _valid_cell_numbering = distributed_data._valid_cell_numbering;
  _valid_edge_numbering = distributed_data._valid_edge_numbering;
  _valid_face_numbering = distributed_data._valid_face_numbering;

  _valid_edge_ownership = distributed_data._valid_edge_ownership;
  _valid_face_ownership = distributed_data._valid_face_ownership;

  _valid_shared_facets_mapping = distributed_data._valid_shared_facets_mapping;

  for (uint i = 0; i < MAX_DIM + 1; i++)
  {
    global_indices[i] = distributed_data.global_indices[i];
    local_indices[i] = distributed_data.local_indices[i];
  }

  for (uint i = 0; i < MAX_DIM; i++)
  {
    adjacent_ranks[i] = distributed_data.adjacent_ranks[i];
    shared[i] = distributed_data.shared[i];
    shared_adj[i] = distributed_data.shared_adj[i];
    ghost[i] = distributed_data.ghost[i];
    ghost_owner[i] = distributed_data.ghost_owner[i];
    shared_mapping[i] = distributed_data.shared_mapping[i];
    ghost_mapping[i] = distributed_data.ghost_mapping[i];
  }

  _num_global_vertex = distributed_data._num_global_vertex;
  _num_global_edge = distributed_data._num_global_edge;
  _num_global_face = distributed_data._num_global_face;
  _num_global_cell = distributed_data._num_global_cell;

  _finalized = distributed_data._finalized;

  _global_vertex_indices_size = distributed_data._global_vertex_indices_size;
  _global_facet_indices_size = distributed_data._global_facet_indices_size;
  _global_cell_indices_size = distributed_data._global_cell_indices_size;

  if (_finalized)
  {
    dolfin_assert(_global_vertex_indices_size > 0);
    dolfin_assert(_global_cell_indices_size > 0);

    _global_vertex_indices = new uint[_global_vertex_indices_size];
    memcpy(_global_vertex_indices, distributed_data._global_vertex_indices,
           _global_vertex_indices_size * sizeof(uint));

    _global_facet_indices = new uint[_global_facet_indices_size];
    memcpy(_global_facet_indices, distributed_data._global_facet_indices,
           _global_facet_indices_size * sizeof(uint));

    _global_cell_indices = new uint[_global_cell_indices_size];
    memcpy(_global_cell_indices, distributed_data._global_cell_indices,
           _global_cell_indices_size * sizeof(uint));
  }

  return *this;
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::empty() const
{
  return (global_indices[0].size() == 0);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::init(uint const dim)
{
  if (dim > 0 && _dim == 0)
  {
    _dim = dim;
    _cell_dim = _dim;
    _facet_dim = _dim - 1;
  }
  else if (dim != _dim)
  {
    error("Trying to set a different topological dimension in distdata.");
  }
}
//-----------------------------------------------------------------------------
void MeshDistributedData::clear()
{
  _dim = 0;
  _cell_dim = 0;
  _facet_dim = 0;
  _max_global_index = 0;

  for (uint i = 0; i < MAX_DIM; ++i)
  {
    adjacent_ranks[i].clear();
    shared[i].clear();
    shared_adj[i].clear();
    ghost[i].clear();
    ghost_owner[i].clear();
    shared_mapping[i].clear();
    ghost_mapping[i].clear();
  }

  for (uint i = 0; i < MAX_DIM + 1; ++i)
  {
    global_indices[i].clear();
    local_indices[i].clear();
  }

  _valid_vertex_numbering = false;
  _valid_edge_numbering = false;
  _valid_face_numbering = false;
  _valid_cell_numbering = false;

  _valid_edge_ownership = false;
  _valid_face_ownership = false;

  _valid_shared_facets_mapping = false;

  delete[] _global_vertex_indices;
  _global_vertex_indices = NULL;

  delete[] _global_facet_indices;
  _global_facet_indices = NULL;

  delete[] _global_cell_indices;
  _global_cell_indices = NULL;

  _finalized = false;

}
//-----------------------------------------------------------------------------
void MeshDistributedData::finalize(uint const dim)
{

  _map<uint, uint>::iterator it;

  if (dim == 0) // Vertices
  {
    delete[] _global_vertex_indices;
    _global_vertex_indices = new uint[global_indices[0].size()];

    for(it = global_indices[0].begin(); it != global_indices[0].end(); ++it)
    {
      _global_vertex_indices[it->first] = it->second;
    }
    _global_vertex_indices_size = global_indices[0].size();
    _max_global_index = _global_vertex_indices_size;
    global_indices[0].clear();
  }
  else if (dim == _facet_dim) // Facets
  {
    delete[] _global_facet_indices;
    _global_facet_indices = new uint[global_indices[dim-1].size()];

    for(it = global_indices[dim-1].begin(); it != global_indices[dim-1].end(); ++it)
    {
      _global_facet_indices[it->first] = it->second;
    }
    _global_facet_indices_size = global_indices[2].size();
    _max_global_index = _global_facet_indices_size;
    global_indices[2].clear();
  }
  else if (dim == _cell_dim) // Cells
  {
    delete[] _global_cell_indices;
    _global_cell_indices = new uint[global_indices[dim].size()];

    for(it = global_indices[dim].begin(); it != global_indices[dim].end(); ++it)
    {
      _global_cell_indices[it->first] = it->second;
    }
    _global_cell_indices_size = global_indices[MAX_DIM].size();
    global_indices[MAX_DIM].clear();
  }
  else
  {
    error("MeshDistributedData::finalize not implemented for %ud.", dim);
  }

  _finalized = true;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_map(uint local_index, uint global_index, uint dim)
{

  if (dim == 0)
  {
    _max_global_index = std::max(_max_global_index, global_index);
  }

  global_indices[dim][local_index] = global_index;
  local_indices[dim][global_index] = local_index;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared(MeshEntity const& m)
{
  set_shared(m.index(), m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared(uint local_index, uint dim)
{
  shared[dim].insert(local_index);
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
  ghost[dim].insert(local_index);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost_owner(MeshEntity const& m, uint rank)
{
  set_ghost_owner(m.index(), rank, m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost_owner(uint i, uint rank, uint dim)
{
  set_shared_adj( i, rank, dim);
  ghost_owner[dim][i] = rank;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared_adj(MeshEntity const& m, uint rank)
{
  set_shared_adj(m.index(), rank, m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared_adj(uint i, uint rank, uint dim)
{
  shared_adj[dim][i].insert(rank);
  adjacent_ranks[dim].insert(rank);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::setall_shared_adj(uint i, _set<uint> const& ranks,
uint dim)
{
  shared_adj[dim][i].clear();
  shared_adj[dim][i].insert(ranks.begin(), ranks.end());
  adjacent_ranks[dim].insert(ranks.begin(), ranks.end());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::setall_shared_adj(MeshEntity const& m,
                                            _set<uint> const& ranks)
{
  setall_shared_adj(m.index(), ranks, m.dim());
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_global(MeshEntity const& e) const
{
  return get_global(e.index(), e.dim());
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_global(uint i, uint dim) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  if (dim == 0 && _finalized)
  {
    return _global_vertex_indices[i];
  }
  else
  {
    dolfin_assert( global_indices[dim].count(i) );
    return global_indices[dim][i];
  }
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_local(MeshEntity const& e) const
{
  return get_local(e.index(), e.dim());
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_local(uint i, uint dim) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  dolfin_assert( local_indices[dim].count(i) );
  return local_indices[dim][i];
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
  }dolfin_assert( ghost_owner[dim].count(local_index) );
  return ghost_owner[dim][local_index];
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
  return shared_adj[dim][local_index];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::num_shared_with(uint rank, uint dim) const
{
  AdjacentMapping::const_iterator it = shared_mapping[dim].find(rank);
  if(it != shared_mapping[dim].end())
  {
    return it->second.first.size();
  }
  return 0;
}

//-----------------------------------------------------------------------------
uint MeshDistributedData::num_ghost_from(uint rank, uint dim) const
{
  AdjacentMapping::const_iterator it = ghost_mapping[dim].find(rank);
  if(it != ghost_mapping[dim].end())
  {
    return it->second.first.size();
  }
  return 0;
}

//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_shared_mapping_to(uint rank,
                                                              uint dim) const
{
  dolfin_assert(shared_mapping[dim].find(rank) != shared_mapping[dim].end());
  return shared_mapping[dim].find(rank)->second.first;
}

//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_shared_mapping_from(uint rank,
                                                                uint dim) const
{
  dolfin_assert(shared_mapping[dim].find(rank) != shared_mapping[dim].end());
  return shared_mapping[dim].find(rank)->second.second;
}

//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_ghost_mapping_to(uint rank,
                                                              uint dim) const
{
  dolfin_assert(shared_mapping[dim].find(rank) != shared_mapping[dim].end());
  return ghost_mapping[dim].find(rank)->second.first;
}

//-----------------------------------------------------------------------------
Array<uint> const& MeshDistributedData::get_ghost_mapping_from(uint rank,
                                                                uint dim) const
{
  dolfin_assert(shared_mapping[dim].find(rank) != shared_mapping[dim].end());
  return ghost_mapping[dim].find(rank)->second.second;
}

//-----------------------------------------------------------------------------
_set<uint> const& MeshDistributedData::get_adj(uint dim) const
{
  dolfin_assert(dim <= _dim);
  return adjacent_ranks[dim];
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
uint MeshDistributedData::get_vertex_global(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  if (_finalized)
  {
    return _global_vertex_indices[i];
  }
  else
  {
    dolfin_assert( global_indices[0].count(i) );
    return global_indices[0][i];
  }

}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_vertex_local(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  dolfin_assert( local_indices[0].count(i) );
  return local_indices[0][i];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_facet_global(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  if (_finalized)
  {
    return _global_facet_indices[i];
  }
  else
  {
    dolfin_assert( _facet_dim != 0 );
    dolfin_assert( global_indices[_facet_dim].count(i) );
    return global_indices[_facet_dim][i];
  }

}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_facet_local(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  dolfin_assert( _facet_dim != 0 );
  dolfin_assert( local_indices[_facet_dim].count(i) );
  return local_indices[_facet_dim][i];
}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_cell_global(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  if (_finalized)
  {
    return _global_cell_indices[i];
  }
  else
  {
    dolfin_assert( _cell_dim != 0 );
    dolfin_assert( global_indices[_cell_dim].count(i) );
    return global_indices[_cell_dim][i];
  }

}
//-----------------------------------------------------------------------------
uint MeshDistributedData::get_cell_local(uint i) const
{
  if (MPI::numProcesses() == 1)
  {
    return i;
  }

  dolfin_assert( _cell_dim != 0 );
  dolfin_assert( local_indices[_cell_dim].count(i) );
  return local_indices[_cell_dim][i];
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::have_global(MeshEntity const& entity) const
{
  return have_global(entity.index(), entity.dim());
}
//-----------------------------------------------------------------------------
bool MeshDistributedData::have_local(MeshEntity const& entity) const
{
  return have_local(entity.index(), entity.dim());
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
  shared_mapping[dim].clear();
  ghost_mapping[dim].clear();
}
//-----------------------------------------------------------------------------
void MeshDistributedData::disp() const
{
  cout << "MeshDistributedData" << endl;
  cout << "-------------------" << endl;

  begin("");
  cout << "Topological dimension     : " << (uint) _dim << endl;
  cout << "Cell dimension            : " << (uint) _cell_dim << endl;
  cout << "Facet dimension           : " << (uint) _facet_dim << endl;
  skip();
  cout << "Maximum global index      : " << (uint) _max_global_index << endl;
  cout << "Number of global vertices : " << (uint) _num_global_vertex << endl;
  cout << "Number of global edges    : " << (uint) _num_global_edge << endl;
  cout << "Number of global faces    : " << (uint) _num_global_face << endl;
  cout << "Number of global cells    : " << (uint) _num_global_cell << endl;
  skip();
  cout << "Valid vertex numbering    : " << (bool) _valid_vertex_numbering
       << endl;
  cout << "Valid edge   numbering    : " << (bool) _valid_edge_numbering
       << endl;
  cout << "Valid face   numbering    : " << (bool) _valid_face_numbering
       << endl;
  cout << "Valid cell   numbering    : " << (bool) _valid_cell_numbering
       << endl;
  skip();
  cout << "Number of shared entities : " << endl;
  for (uint d = 0; d < _dim; ++d)
  {
    cout << "  - dim " << " : " << (uint) this->num_shared(d) << endl;
  }
  skip();
  cout << "Number of ghost entities : " << endl;
  for (uint d = 0; d < _dim; ++d)
  {
    cout << "  - dim " << " : " << (uint) this->num_ghost(d) << endl;
  }
  end();
}

}

