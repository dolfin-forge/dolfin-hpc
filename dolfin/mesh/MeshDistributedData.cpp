// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#include "MeshDistributedData.h"
#include "MeshEntity.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include <dolfin/main/MPI.h>
#include <dolfin/log/log.h>
#include <mpi.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
MeshDistributedData::MeshDistributedData() : _size(0), _cell_size(0),
					     _max_global_index(0),
					     _num_shared(0), _num_ghost(0), 
					     _valid_vertex_numbering(false),
					     _valid_cell_numbering(false),
					     _valid_edge_numbering(false),
					     _valid_face_numbering(false)
					    					     
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

  _size = distributed_data._size;
  _cell_size = distributed_data._cell_size;
  _num_shared = distributed_data._num_shared;
  _num_ghost = distributed_data._num_ghost;

  _valid_vertex_numbering = distributed_data._valid_vertex_numbering;
  _valid_cell_numbering = distributed_data._valid_cell_numbering;
  _valid_edge_numbering = distributed_data._valid_edge_numbering;
  _valid_face_numbering = distributed_data._valid_face_numbering;

  _max_global_index = distributed_data._max_global_index;
  
  global_vertex_indices = distributed_data.global_vertex_indices;
  local_vertex_indices = distributed_data.local_vertex_indices;

  global_edge_indices = distributed_data.global_edge_indices;
  local_edge_indices = distributed_data.local_edge_indices;

  global_face_indices = distributed_data.global_face_indices;
  local_face_indices = distributed_data.local_face_indices;

  global_cell_indices = distributed_data.global_cell_indices;
  local_cell_indices = distributed_data.local_cell_indices;

  shared_vertices = distributed_data.shared_vertices;
  ghost_vertices = distributed_data.ghost_vertices;
  ghost_owner = distributed_data.ghost_owner;

  return *this;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::clear()
{
  _size = _cell_size = _num_shared = _num_ghost = _max_global_index = 0;
  shared_vertices.clear(); 
  ghost_vertices.clear();
  ghost_owner.clear();

  global_vertex_indices.clear();
  local_vertex_indices.clear();

  global_edge_indices.clear();
  local_edge_indices.clear();

  global_face_indices.clear();
  local_face_indices.clear();

  global_cell_indices.clear();
  local_cell_indices.clear();

  _valid_vertex_numbering = _valid_cell_numbering = false;
  _valid_edge_numbering = _valid_face_numbering = false;

}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_map(uint local_index, uint global_index, uint dim)
{

  switch(dim) {
  case 0:
    if(global_vertex_indices.count(local_index) == 0) {
      _size++;
      _max_global_index = std::max(_max_global_index, global_index);
    }
    global_vertex_indices[ local_index ] = global_index;
    local_vertex_indices[ global_index ] = local_index;
    break;
  case 1:
    global_edge_indices[ local_index ] = global_index;
    local_edge_indices[ global_index ] = local_index;
    break;
  case 2:
    global_face_indices[ local_index ] = global_index;
    local_face_indices[ global_index ] = local_index;
    break;
  }      
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared(Vertex& v)
{
  set_shared(v.index());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_shared(uint local_index)
{
  if( shared_vertices.count(local_index) == 0) {
    _num_shared++;
    shared_vertices.insert(local_index);
  }
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost(Vertex& v)
{
  set_ghost(v.index());
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost(uint local_index)
{
  if(ghost_vertices.count(local_index) == 0) {
    _num_ghost++;
    ghost_vertices.insert(local_index);
    set_shared(local_index);
  }
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost_owner(Vertex& v, uint rank)
{
  set_ghost_owner(v.index(), rank);
}
//-----------------------------------------------------------------------------
void MeshDistributedData::set_ghost_owner(uint i, uint rank)
{
  ghost_owner[i] = rank;
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

  switch(dim) {
  case 0:    
    dolfin_assert( global_vertex_indices.count(i) );
    return global_vertex_indices[i];
    break;
  case 1:
    dolfin_assert(global_edge_indices.count(i) );
    return global_edge_indices[i];
    break;
  case 2:
    dolfin_assert( global_face_indices.count(i) );
    return global_face_indices[i];
    break;    
  }
  
  return 0;

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

  switch(dim) {
  case 0:
    dolfin_assert( local_vertex_indices.count(i) );
    return local_vertex_indices[i];
    break;
  case 1:
    dolfin_assert( local_edge_indices.count(i) );
    return local_edge_indices[i];
    break;
  case 2:
    dolfin_assert( local_face_indices.count(i) );
    return local_face_indices[i];
    break;
  }

  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_owner(Vertex& v) 
{
  return get_owner(v.index());
}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_owner(uint local_index) 
{ 
  dolfin_assert( ghost_owner.count(local_index) );
  return ghost_owner[local_index];
}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_cell_global(uint i)
{
  if(MPI::numProcesses() == 1) 
    return i;
  dolfin_assert( global_cell_indices.count(i) );
  return global_cell_indices[i];

}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_cell_local(uint i)
{
  if(MPI::numProcesses() == 1) 
    return i;
  dolfin_assert( local_cell_indices.count(i) );
  return local_cell_indices[i];
}
//-----------------------------------------------------------------------------
void MeshDistributedData::invalid_numbering() 
{
  _valid_vertex_numbering = _valid_cell_numbering = false;
  _valid_edge_numbering = _valid_face_numbering = false;
}
//-----------------------------------------------------------------------------


