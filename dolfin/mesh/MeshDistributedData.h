// Copyright (C) 2008 Niclas Jansson. 
// Licensed under the GNU LGPL Version 2.1. 
//


#ifndef __MESH_DISTRIBUTED_DATA_H
#define __MESH_DISTRIBUTED_DATA_H

#include <dolfin/common/types.h>
#include <dolfin/main/MPI.h>
#include <dolfin/log/log.h>
#include <set>
#include <map>

namespace dolfin
{
  class Vertex;
  class Edge;
  class Face;
  class MeshEntity;
  class MeshDistributedData
  {
  public:
    MeshDistributedData();
    
    ~MeshDistributedData();

    const MeshDistributedData& operator=(const MeshDistributedData& distributed_data);
    
    void clear();
    void set_map(uint local_index, uint global_index, uint dim);
    
    void set_shared(uint local_index);
    void set_ghost(uint local_index);

    void set_shared(Vertex& v);
    void set_ghost(Vertex& v);
    
    uint get_global(uint i, uint dim);
    uint get_global(MeshEntity& e);

    uint get_local(uint i, uint dim);
    uint get_local(MeshEntity& e);
    
    uint get_cell_global(uint i);
    uint get_cell_local(uint i);

    void set_ghost_owner(uint i, uint rank);
    void set_ghost_owner(Vertex& v, uint rank);

    inline void set_global_numVertices(uint num_global) 
    { _num_global_vertex = num_global; }

    inline void set_global_numEdges(uint num_global) 
    { _num_global_edge = num_global; }

    inline void set_global_numFaces(uint num_global) 
    { _num_global_face = num_global; }

    inline void set_global_numCells(uint num_global) 
    { _num_global_cell = num_global; }
    
    void invalid_numbering();

    uint get_owner(uint local_index);
    uint get_owner(Vertex& v);

    inline bool have_global(uint i) 
    {return (MPI::numProcesses() > 1 ? (local_vertex_indices.count(i) > 0) : true);}
    
    inline bool have_local(uint i) 
    {return (MPI::numProcesses() > 1 ? (global_vertex_indices.count(i) > 0) : true);}
    
    inline bool is_shared(uint i)
    {return (MPI::numProcesses() > 1 ? (shared_vertices.count(i) > 0) : true);}

    inline bool is_ghost(uint i)
    {return (MPI::numProcesses() > 1 ? (ghost_vertices.count(i) > 0) : true);}

    inline uint num_shared() { return _num_shared; }
    inline uint num_ghost() { return _num_ghost; }

    inline uint global_numVertices() { return _num_global_vertex; }
    
    inline uint global_numEdges() { return _num_global_edge; }

    inline uint global_numFaces() { return _num_global_face; }

    inline uint global_numCells() { return _num_global_cell; }
    
    inline uint max_index() { return _max_global_index; }


  private:
    
    uint _size, _cell_size,_max_global_index;
    uint _num_shared, _num_ghost;
    uint _num_global_vertex, _num_global_edge;
    uint _num_global_face, _num_global_cell;


    bool _valid_vertex_numbering, _valid_cell_numbering,
      _valid_edge_numbering, _valid_face_numbering;

    std::map<uint, uint> global_vertex_indices;
    std::map<uint, uint> local_vertex_indices;

    std::map<uint, uint> global_edge_indices;
    std::map<uint, uint> local_edge_indices;

    std::map<uint, uint> global_face_indices;
    std::map<uint, uint> local_face_indices;

    std::map<uint, uint> global_cell_indices;
    std::map<uint, uint> local_cell_indices;
  
    std::map<uint, uint> ghost_owner;
    
    std::set<uint> shared_vertices;
    std::set<uint> ghost_vertices;

    friend class MeshGhostIterator;
    friend class MeshSharedIterator;
    friend class MeshRenumber;

   };
  
  class MeshGhostIterator 
  {
  public:
  MeshGhostIterator(MeshDistributedData& distdata) : _distdata(distdata) 
    { _iter = _distdata.ghost_vertices.begin();}
    
    ~MeshGhostIterator() {}
    MeshGhostIterator& operator++() { ++_iter; return *this;}
    inline uint index() const { return *_iter; }
    inline uint owner() { return _distdata.get_owner(*_iter); }   
    inline bool end() const { return _iter == _distdata.ghost_vertices.end();}
  private:
    MeshDistributedData& _distdata;
    std::set<uint>::iterator _iter;
  };

  class MeshSharedIterator 
  {
  public:
  MeshSharedIterator(MeshDistributedData& distdata) : _distdata(distdata) 
    { _iter = _distdata.shared_vertices.begin(); }
      
    ~MeshSharedIterator() {}
    MeshSharedIterator& operator++() { ++_iter; return *this;}
    inline uint index() const { return *_iter; }
    inline bool end() const { return _iter == _distdata.shared_vertices.end();}
   
  private:
    MeshDistributedData& _distdata;
    std::set<uint>::iterator _iter;
  };
  
}

#endif
