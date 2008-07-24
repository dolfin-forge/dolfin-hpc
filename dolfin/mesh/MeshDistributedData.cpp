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
MeshDistributedData::MeshDistributedData() : _max_global_index(0),
					     _valid_vertex_numbering(false),
					     _valid_cell_numbering(false),
					     _valid_edge_numbering(false),
					     _valid_face_numbering(false),
					     _valid_edge_ownership(false),
					     _valid_face_ownership(false)
					    					     
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
    ghost[i] = distributed_data.ghost[i];
    ghost_owner[i] = distributed_data.ghost_owner[i];
  }

  _num_global_vertex = distributed_data._num_global_vertex;
  _num_global_edge = distributed_data._num_global_edge;
  _num_global_face = distributed_data._num_global_face;
  _num_global_cell = distributed_data._num_global_cell;


  return *this;
}
//-----------------------------------------------------------------------------
void MeshDistributedData::clear()
{
  _max_global_index = 0;

  for(uint i = 0; i < 3; i++) {
    shared[i].clear(); 
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
dolfin::uint MeshDistributedData::get_global(MeshEntity& e)
{
  return get_global( e.index(), e.dim());
}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_global(uint i, uint dim)
{
  if(MPI::numProcesses() == 1) 
    return i;

  dolfin_assert( global_indices[dim].count(i) );
  return global_indices[dim][i];

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
  dolfin_assert( ghost_owner[dim].count(local_index) );
  return ghost_owner[dim][local_index];
}
//-----------------------------------------------------------------------------
dolfin::uint MeshDistributedData::get_cell_global(uint i)
{
  if(MPI::numProcesses() == 1) 
    return i;
  dolfin_assert( global_indices[3].count(i) );
  return global_indices[3][i];

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


