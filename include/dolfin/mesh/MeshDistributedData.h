// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __MESH_DISTRIBUTED_DATA_H
#define __MESH_DISTRIBUTED_DATA_H

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

class Vertex;
class Edge;
class Face;
class MeshEntity;

class MeshDistributedData
{

  // Friends
  friend class MeshGhostIterator;
  friend class MeshSharedIterator;
  friend class MeshRenumber;

public:

  MeshDistributedData();

  ~MeshDistributedData();

  const MeshDistributedData& operator=(
      const MeshDistributedData& distributed_data);

  void finalize(uint dim);

  void set_map(uint local_index, uint global_index, uint dim);

  void set_shared(uint local_index, uint dim);
  void set_ghost(uint local_index, uint dim);

  void set_shared(MeshEntity& m);
  void set_ghost(MeshEntity& m);

  uint get_global(uint i, uint dim);
  uint get_global(MeshEntity& e);

  uint get_local(uint i, uint dim);
  uint get_local(MeshEntity& e);

  uint get_cell_global(uint i);
  uint get_cell_local(uint i);

  void set_ghost_owner(uint i, uint rank, uint dim);
  void set_ghost_owner(MeshEntity& m, uint rank);

  void set_shared_adj(uint i, uint rank, uint dim);
  void set_shared_adj(MeshEntity& m, uint rank);

  _set<uint>& get_shared_adj(uint local_index, uint dim);
  _set<uint>& get_shared_adj(MeshEntity& m);

  void set_global_numVertices(uint num_global);

  void set_global_numEdges(uint num_global);

  void set_global_numFaces(uint num_global);

  void set_global_numCells(uint num_global);

  void invalid_numbering();

  void invalid_ownership();

  uint get_owner(uint local_index, uint dim);
  uint get_owner(MeshEntity& m);
  void remap_owner(int* mapping);

  bool have_global(uint i, uint dim);

  bool have_local(uint i, uint dim);

  bool is_shared(uint i, uint dim);

  bool is_ghost(uint i, uint dim);

  uint num_shared(uint dim);

  uint num_ghost(uint dim);

  uint global_numVertices();

  uint global_numEdges();

  uint global_numFaces();

  uint global_numCells();

  uint max_index();

  void flush_edges();

  void flush_faces();

private:

  void clear();

  //--- ATTRIBUTES ------------------------------------------------------------

  uint _max_global_index;
  uint _num_global_vertex;
  uint _num_global_edge;
  uint _num_global_face;
  uint _num_global_cell;

  bool _valid_vertex_numbering;
  bool _valid_cell_numbering;
  bool _valid_edge_numbering;
  bool _valid_face_numbering;

  bool _valid_edge_ownership;
  bool _valid_face_ownership;

  _map<uint, uint> global_indices[4];
  _map<uint, uint> local_indices[4];

  _map<uint, uint> ghost_owner[3];
  _map<uint, _set<uint> > shared_adj[3];

  _set<uint> shared[3];
  _set<uint> ghost[3];

  bool finalized;
  uint *_global_indices;
  uint *_global_cell_indices;

  uint _global_indices_size;
  uint _global_cell_indices_size;

};

class MeshGhostIterator
{

public:

  MeshGhostIterator(MeshDistributedData& distdata, uint i);

  ~MeshGhostIterator();

  MeshGhostIterator& operator++();

  uint index() const;

  uint owner();

  bool end() const;

private:

  MeshDistributedData& _distdata;

  _set<uint>::iterator _iter;
  uint _dim;
};

class MeshSharedIterator
{

public:

  MeshSharedIterator(MeshDistributedData& distdata, uint i);

  ~MeshSharedIterator();

  MeshSharedIterator& operator++();

  uint index() const;

  bool end() const;

  _set<uint> adj() const;

private:

  MeshDistributedData& _distdata;

  _set<uint>::iterator _iter;
  uint _dim;
};

//--- INLINES -----------------------------------------------------------------

//-----------------------------------------------------------------------------
inline void MeshDistributedData::set_global_numVertices(uint num_global)
{
  _num_global_vertex = num_global;
}

//-----------------------------------------------------------------------------
inline void MeshDistributedData::set_global_numEdges(uint num_global)
{
  _num_global_edge = num_global;
}

//-----------------------------------------------------------------------------
inline void MeshDistributedData::set_global_numFaces(uint num_global)
{
  _num_global_face = num_global;
}

//-----------------------------------------------------------------------------
inline void MeshDistributedData::set_global_numCells(uint num_global)
{
  _num_global_cell = num_global;
}

//-----------------------------------------------------------------------------
inline void MeshDistributedData::invalid_numbering()
{
  _valid_vertex_numbering = false;
  _valid_cell_numbering = false;
  _valid_edge_numbering = false;
  _valid_face_numbering = false;
  finalized = false;
}

//-----------------------------------------------------------------------------
inline void MeshDistributedData::invalid_ownership()
{
  _valid_edge_numbering = false;
  _valid_face_numbering = false;
  flush_edges();
  flush_faces();
  finalized = false;
}

//-----------------------------------------------------------------------------
inline bool MeshDistributedData::have_global(uint i, uint dim)
{
  return (MPI::numProcesses() > 1 ? (local_indices[dim].count(i) > 0) : true);
}

//-----------------------------------------------------------------------------
inline bool MeshDistributedData::have_local(uint i, uint dim)
{
  return (MPI::numProcesses() > 1 ? (global_indices[dim].count(i) > 0) : true);
}

//-----------------------------------------------------------------------------
inline bool MeshDistributedData::is_shared(uint i, uint dim)
{
  return (MPI::numProcesses() > 1 ? (shared[dim].count(i) > 0) : true);
}

//-----------------------------------------------------------------------------
inline bool MeshDistributedData::is_ghost(uint i, uint dim)
{
  return (MPI::numProcesses() > 1 ? (ghost[dim].count(i) > 0) : false);
}

//-----------------------------------------------------------------------------
inline uint MeshDistributedData::num_shared(uint dim)
{
  return shared[dim].size();
}

//-----------------------------------------------------------------------------
inline uint MeshDistributedData::num_ghost(uint dim)
{
  return ghost[dim].size();
}

//-----------------------------------------------------------------------------
inline uint MeshDistributedData::global_numVertices()
{
  return _num_global_vertex;
}

//-----------------------------------------------------------------------------
inline uint MeshDistributedData::global_numEdges()
{
  return _num_global_edge;
}

//-----------------------------------------------------------------------------
inline uint MeshDistributedData::global_numFaces()
{
  return _num_global_face;
}

//-----------------------------------------------------------------------------
inline uint MeshDistributedData::global_numCells()
{
  return _num_global_cell;
}

//-----------------------------------------------------------------------------
inline uint MeshDistributedData::max_index()
{
  return _max_global_index;
}

//-----------------------------------------------------------------------------
inline void MeshDistributedData::flush_edges()
{
  shared[1].clear();
  ghost[1].clear();
  ghost_owner[1].clear();
}

//-----------------------------------------------------------------------------
inline void MeshDistributedData::flush_faces()
{
  shared[2].clear();
  ghost[2].clear();
  ghost_owner[2].clear();
}

//-----------------------------------------------------------------------------
inline uint MeshGhostIterator::index() const
{
  return *_iter;
}

//-----------------------------------------------------------------------------
inline uint MeshGhostIterator::owner()
{
  return _distdata.get_owner(*_iter, _dim);
}

//-----------------------------------------------------------------------------
inline bool MeshGhostIterator::end() const
{
  return _iter == _distdata.ghost[_dim].end();
}

//-----------------------------------------------------------------------------
inline uint MeshSharedIterator::index() const
{
  return *_iter;
}

//-----------------------------------------------------------------------------
inline bool MeshSharedIterator::end() const
{
  return _iter == _distdata.shared[_dim].end();
}

//-----------------------------------------------------------------------------
inline _set<uint> MeshSharedIterator::adj() const
{
  return _distdata.shared_adj[_dim][*_iter];
}

//-----------------------------------------------------------------------------

}

#endif
