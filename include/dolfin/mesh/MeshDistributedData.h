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
class Mesh;
class MeshEntity;
class MeshDistributedData
{
  static uint const MAX_DIM = 3;

public:

  MeshDistributedData(Mesh const& mesh);

  ~MeshDistributedData();

  const MeshDistributedData& operator=(const MeshDistributedData& distributed_data);

  void clear();
  void finalize(uint dim);

  //--- These cannot be const as they access the element in the map with [].

  uint get_global(uint i, uint dim);
  uint get_global(MeshEntity& e);

  uint get_local(uint i, uint dim);
  uint get_local(MeshEntity& e);

  uint get_vertex_global(uint i);
  uint get_vertex_local(uint i);

  uint get_facet_global(uint i);
  uint get_facet_local(uint i);

  uint get_cell_global(uint i);
  uint get_cell_local(uint i);

  //---

  void set_map(uint local_index, uint global_index, uint dim);

  void set_shared(uint local_index, uint dim);
  void set_ghost(uint local_index, uint dim);

  void set_shared(MeshEntity& m);
  void set_ghost(MeshEntity& m);

  void set_ghost_owner(uint i, uint rank, uint dim);
  void set_ghost_owner(MeshEntity& m, uint rank);

  void set_shared_adj(uint i, uint rank, uint dim);
  void set_shared_adj(MeshEntity& m, uint rank);

  _set<uint>& get_shared_adj(uint local_index, uint dim);
  _set<uint>& get_shared_adj(MeshEntity& m);

  inline void set_global_numVertices(uint num_global)
  { _num_global_vertex = num_global; }

  inline void set_global_numEdges(uint num_global)
  { _num_global_edge = num_global; }

  inline void set_global_numFaces(uint num_global)
  { _num_global_face = num_global; }

  inline void set_global_numCells(uint num_global)
  { _num_global_cell = num_global; }

  inline void invalid_numbering()
  {
    _valid_vertex_numbering = _valid_cell_numbering = false;
    _valid_edge_numbering = _valid_face_numbering = false;
    finalized = false;
  }

  inline void invalid_ownership()
  {
    _valid_edge_numbering = _valid_face_numbering = false;
    flush_edges(); flush_faces();
    finalized = false;
  }

  uint get_owner(uint local_index, uint dim);
  uint get_owner(MeshEntity& m);
  void remap_owner(int* mapping);

  inline bool have_global(uint i, uint dim) const
  {return (MPI::numProcesses() > 1 ? (local_indices[dim].count(i) > 0) : true);}

  inline bool have_local(uint i, uint dim) const
  {return (MPI::numProcesses() > 1 ? (global_indices[dim].count(i) > 0) : true);}

  inline bool is_shared(uint i, uint dim) const
  {return (MPI::numProcesses() > 1 ? (shared[dim].count(i) > 0) : true);}

  inline bool is_ghost(uint i, uint dim) const
  {return (MPI::numProcesses() > 1 ? (ghost[dim].count(i) > 0) : false);}

  inline uint num_shared(uint dim) const {return shared[dim].size(); }

  inline uint num_ghost(uint dim) const {return ghost[dim].size(); }

  inline uint global_numVertices() const { return _num_global_vertex; }

  inline uint global_numEdges() const { return _num_global_edge; }

  inline uint global_numFaces() const { return _num_global_face; }

  inline uint global_numCells() const { return _num_global_cell; }

  inline uint max_index() const { return _max_global_index; }

  inline void flush_edges()
  { shared[1].clear(); ghost[1].clear(); ghost_owner[1].clear();}

  inline void flush_faces()
  { shared[2].clear(); ghost[2].clear(); ghost_owner[2].clear();}

protected:

private:

  Mesh const& _mesh;

  uint _max_global_index;
  uint _num_global_vertex;
  uint _num_global_edge;
  uint _num_global_face;
  uint _num_global_cell;

  bool _valid_vertex_numbering;
  bool _valid_edge_numbering;
  bool _valid_face_numbering;
  bool _valid_cell_numbering;

  bool _valid_edge_ownership;
  bool _valid_face_ownership;

  _map<uint, uint> global_indices[MAX_DIM+1];
  _map<uint, uint> local_indices[MAX_DIM+1];

  _map<uint, uint> ghost_owner[MAX_DIM];
  _map<uint, _set<uint> > shared_adj[MAX_DIM];

  _set<uint> shared[MAX_DIM];
  _set<uint> ghost[MAX_DIM];

  bool finalized;
  uint *_global_vertex_indices;
  uint *_global_facet_indices;
  uint *_global_cell_indices;

  uint _global_vertex_indices_size;
  uint _global_facet_indices_size;
  uint _global_cell_indices_size;

  friend class MeshGhostIterator;
  friend class MeshSharedIterator;
  friend class MeshRenumber;

 };

class MeshGhostIterator
{
public:
  MeshGhostIterator(MeshDistributedData& distdata, uint i) : _distdata(distdata)
  { _iter = _distdata.ghost[i].begin(); _dim = i;}

  ~MeshGhostIterator() {}
  MeshGhostIterator& operator++() { ++_iter; return *this;}
//    MeshGhostIterator operator++(int) { _iter++; return *this;}
  inline uint index() const { return *_iter; }
  inline uint owner() { return _distdata.get_owner(*_iter, _dim); }
  inline bool end() const { return _iter == _distdata.ghost[_dim].end();}

private:
  MeshDistributedData& _distdata;
  _set<uint>::iterator _iter;
  uint _dim;
};

class MeshSharedIterator
{
public:
  MeshSharedIterator(MeshDistributedData& distdata, uint i) : _distdata(distdata)
  { _iter = _distdata.shared[i].begin(); _dim = i;}

  ~MeshSharedIterator() {}
  MeshSharedIterator& operator++() { ++_iter; return *this;}
//    MeshSharedIterator operator++(int) { _iter++; return *this;}
  inline uint index() const { return *_iter; }
  inline bool end() const { return _iter == _distdata.shared[_dim].end();}
  inline _set<uint> adj() const { return _distdata.shared_adj[_dim][*_iter]; }

private:
  MeshDistributedData& _distdata;
  _set<uint>::iterator _iter;
  uint _dim;
};

}

#endif
