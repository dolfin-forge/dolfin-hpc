// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/config/dolfin_config.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Face.h>
#include <dolfin/main/MPI.h>
#include <dolfin/log/log.h>
#include <string.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
MeshDistributedData::MeshDistributedData() : _max_global_index(0),
					     _valid_vertex_numbering(false),
					     _valid_cell_numbering(false),
					     _valid_edge_numbering(false),
					     _valid_face_numbering(false),
					     _valid_edge_ownership(false),
					     _valid_face_ownership(false),
					     finalized(false),
					     _global_indices(0),
					     _global_cell_indices(0)
					    					     
{

}
//-----------------------------------------------------------------------------
MeshDistributedData::~MeshDistributedData()
{ 
  clear();
}
//-----------------------------------------------------------------------------
const MeshDistributedData& MeshDistributedData::operator=(const MeshDistributedData& distributed_data)
{
  clear();

  _max_global_index = distributed_data._max_global_index;

  _valid_vertex_numbering = distributed_data._valid_vertex_numbering;
  _valid_cell_numbering = distributed_data._valid_cell_numbering;
  _valid_edge_numbering = distributed_data._valid_edge_numbering;
  _valid_face_numbering = distributed_data._valid_face_numbering;

  _valid_edge_ownership = distributed_data._valid_edge_ownership;
  _valid_face_ownership = distributed_data._valid_face_ownership;


  for(uint i = 0 ; i < 4; i++) {
    global_indices[i] = distributed_data.global_indices[i];
    local_indices[i] = distributed_data.local_indices[i];
  }
 
  for(uint i = 0 ; i < 3; i++) {
    shared[i] = distributed_data.shared[i];
    shared_adj[i] = distributed_data.shared_adj[i];
    ghost[i] = distributed_data.ghost[i];
    ghost_owner[i] = distributed_data.ghost_owner[i];
  }

  _num_global_vertex = distributed_data._num_global_vertex;
  _num_global_edge = distributed_data._num_global_edge;
  _num_global_face = distributed_data._num_global_face;
  _num_global_cell = distributed_data._num_global_cell;

  finalized =  distributed_data.finalized;

  _global_indices_size = distributed_data._global_indices_size;
  _global_cell_indices_size = distributed_data._global_cell_indices_size;

  if ( finalized ) {
    dolfin_assert(_global_indices_size > 0);
    dolfin_assert(_global_cell_indices_size > 0);

    _global_indices = new uint[_global_indices_size];
    memcpy(_global_indices, distributed_data._global_indices, 
	   _global_indices_size * sizeof(uint));

    _global_cell_indices = new uint[_global_cell_indices_size];
    memcpy(_global_cell_indices, distributed_data._global_cell_indices, 
	   _global_cell_indices_size * sizeof(uint));
  }

    
  return *this;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::clear()
{
  _max_global_index = 0;

  for(uint i = 0; i < 3; i++) {
    shared[i].clear(); 
    shared_adj[i].clear();
    ghost[i].clear();
    ghost_owner[i].clear();
  }

  for(uint i = 0; i < 4; i++) {
    global_indices[i].clear();
    local_indices[i].clear();
  }

  _valid_vertex_numbering = _valid_cell_numbering = false;
  _valid_edge_numbering = _valid_face_numbering = false;

  _valid_edge_ownership = _valid_face_ownership = false;

  if( _global_indices )
    delete[] _global_indices;
  _global_indices = 0;

  if( _global_cell_indices )
    delete[] _global_cell_indices;
  _global_cell_indices = 0;

  finalized = false;

}
//-----------------------------------------------------------------------------
void MeshDistributedData::finalize(uint dim)
{

  _map<uint, uint>::iterator it;

  switch(dim)
  {
  case 0:
    if(_global_indices)
      delete[] _global_indices;    
    _global_indices = new uint[global_indices[0].size()];

    for(it = global_indices[0].begin(); it != global_indices[0].end(); ++it) 
      _global_indices[it->first] = it->second;
    _global_indices_size = global_indices[0].size();
    _max_global_index = _global_indices_size;
    global_indices[0].clear();
    break;
  case 3:
    if(_global_cell_indices)
      delete[] _global_cell_indices;
    
    _global_cell_indices = new uint[global_indices[dim].size()];
    for(it = global_indices[dim].begin(); it != global_indices[dim].end(); ++it)
      _global_cell_indices[it->first] = it->second;
    _global_cell_indices_size = global_indices[3].size();
    global_indices[3].clear();
    break;
  default:
    warning("Not implemented yet!");
    break;
  }
  
  finalized = true;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_map(uint local_index, uint global_index, uint dim)
{

  if( dim == 0) 
    _max_global_index = std::max(_max_global_index, global_index);

  global_indices[dim][ local_index ] = global_index;
  local_indices[dim][ global_index ] = local_index;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared(MeshEntity& m)
{
  set_shared(m.index(), m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared(uint local_index, uint dim)
{
  shared[dim].insert(local_index);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost(MeshEntity& m)
{
  set_ghost(m.index(), m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost(uint local_index, uint dim)
{
  ghost[dim].insert(local_index);
  set_shared(local_index, dim);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost_owner(MeshEntity& m, uint rank)
{
  set_ghost_owner(m.index(), rank, m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost_owner(uint i, uint rank, uint dim)
{
  ghost_owner[dim][i] = rank;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared_adj(MeshEntity& m, uint rank)
{
  set_shared_adj(m.index(), rank, m.dim());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared_adj(uint i, uint rank, uint dim)
{
  shared_adj[dim][i].insert(rank);
}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_global(MeshEntity& e)
{
  return get_global( e.index(), e.dim());
}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_global(uint i, uint dim)
{
  if(MPI::numProcesses() == 1) 
    return i;

  if( dim == 0 && finalized) {
    return _global_indices[i];
  }
  else {
    dolfin_assert( global_indices[dim].count(i) );  
    return global_indices[dim][i];
  }
}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_local(MeshEntity& e)
{
  return get_local(e.index(), e.dim());
}
//-----------------------------------------------------------------------------
 dolfin::uint MeshDistributedData::get_local(uint i, uint dim)
{
  if(MPI::numProcesses() == 1)
    return i;

  dolfin_assert( local_indices[dim].count(i) );
  return local_indices[dim][i];
}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_owner(MeshEntity& e) 
{
  return get_owner(e.index(), e.dim());
}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_owner(uint local_index, uint dim) 
{ 
  if(MPI::numProcesses() == 1)
    return 0;
  dolfin_assert( ghost_owner[dim].count(local_index) );
  return ghost_owner[dim][local_index];
}
//-----------------------------------------------------------------------------
_set<dolfin::uint>& MeshDistributedData::get_shared_adj(MeshEntity& m)
{
  return get_shared_adj(m.index(), m.dim());
}
//-----------------------------------------------------------------------------
_set<dolfin::uint>& MeshDistributedData::get_shared_adj(uint local_index, 
						       uint dim) 
{
  dolfin_assert(is_shared(local_index, dim));
  dolfin_assert(!is_ghost(local_index, dim));

  return shared_adj[dim][local_index];
}
//-----------------------------------------------------------------------------
void MeshDistributedData::remap_owner(int* mapping) 
{

  for (uint i = 0;  i < 3; i++)
  {
    for (MeshGhostIterator it(*this, i); !it.end(); ++it)
      set_ghost_owner(it.index(), mapping[it.owner()], i);
#ifdef ENABLE_P1_OPTIMIZATIONS    
    break;
#endif
  }
  
}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_cell_global(uint i)
{
  if(MPI::numProcesses() == 1) 
    return i;

  if ( finalized )
    return _global_cell_indices[i];
  else {
    dolfin_assert( global_indices[3].count(i) );
    return global_indices[3][i];
  }

}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_cell_local(uint i)
{
  if(MPI::numProcesses() == 1) 
    return i;
  dolfin_assert( local_indices[3].count(i) );
  return local_indices[3][i];
}
//-----------------------------------------------------------------------------


