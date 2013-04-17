// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2010.
// Modified by Balthasar Reuter, 2013
//

#ifndef __D_MESH_H
#define __D_MESH_H

#include <dolfin/common/types.h>
#include <dolfin/main/MPI.h>

#include <vector>
#include <list>
#include <set>
#include <map>

namespace dolfin
{
  class Mesh;
  class DCell;
  class DVertex;
  class CellType;

  class DMesh
  {
  public:
    DMesh();
    ~DMesh();

    typedef std::pair<int, int> EdgeKey;    

    typedef struct __edge__ {
      uint mv;
      uint v1; 
      uint v2;
      uint owner;
    } prop_edge;

    typedef std::pair<uint, prop_edge> Propagation;

    void addVertex(DVertex* v);
    
    void addCell(DCell* c, std::vector<DVertex*> vs, int parent_id);
    void removeCell(DCell* c);

    void imp(Mesh& mesh);
    void exp(Mesh& mesh);
    void number();

    void bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1);

    void bisectMarked(std::vector<bool> marked_ids);

    DCell* opposite(DCell* dcell, DVertex* v1, DVertex* v2);

    inline void propagate_refinement(std::vector<Propagation>& propagated,
                                     bool& empty)
    {
      if ( MPI::numProcesses() & ( MPI::numProcesses() - 1 ) )
        propagate_naive(propagated, empty); 
      else
        propagate_hypercube(propagated, empty); 
    }

    void propagate_naive(std::vector<Propagation>& propagated, bool& empty);

    void propagate_hypercube(std::vector<Propagation>& propagated, bool& empty);
    
    std::set<DVertex *> vertices;
    std::list<DCell *> cells;

    CellType* cell_type;
    uint d;
    uint glb_max;
    uint _salt;

    /// Start offset for new global id
    uint _start_offset, _max;

    /// Propagation buffer
    std::vector<Propagation> propagate;

    /// Map between global number of boundary vertex to vertex
    _map<uint, DVertex*> bc_dvs;
#if (__GNUG__ || __sgi )
    std::map<EdgeKey, DVertex*> ref_edge;
#else
    _map<EdgeKey, DVertex*> ref_edge;
#endif

    /// Comparison operator for index/value pairs
    struct less_pair : public std::binary_function<std::pair<uint, prop_edge>,
                       std::pair<uint, prop_edge>, bool>
    {
      bool operator()(std::pair<uint, prop_edge> x, std::pair<uint, prop_edge> y)
      {
        return x.first < y.first;
      }
    };


    /// Construct a edge id from given vertices
    inline EdgeKey edge_key(int id1, int id2) 
    {
      dolfin_assert( id2 != id1 );
      if(id2 < id1)
      {
        EdgeKey key(id2,id1);    
        return key;
      }
      else 
      {
        EdgeKey key(id1,id2);    
        return key;
      }      
    }

  };
}

#endif