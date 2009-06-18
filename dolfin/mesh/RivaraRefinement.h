// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009.
//

#ifndef __RIVARAREFINEMENT_H
#define __RIVARAREFINEMENT_H

#include <list>
#include <vector>

#include <dolfin/mesh/MeshFunction.h>

namespace dolfin
{

  class Cell;
  class Edge;
  class Mesh;

  class RivaraRefinement
  {
  public:
    
    /// Refine simplicial mesh locally by recursive edge bisection 
    static void refine(Mesh& mesh, 
		       MeshFunction<bool>& cell_marker,
		       real tf = 0.0, 
		       real tb = 0.0, 
		       real ts = 0.0);
  };

  class DVertex;
  class DCell;
  class DMesh;
  class DEdge;
  
  class DVertex
  {
  public:
    DVertex();

    int id;
    int glb_id;

    std::list<DCell *> cells;
    Point p;

    bool on_boundary;
    bool shared;
    bool ghosted;

    int owner;
  };
    
  class DCell
  {
  public:
    DCell();
    bool has_edge(DVertex* v1, DVertex* v2);

    int id;
    int parent_id;

    std::vector<DVertex *> vertices;

    bool deleted;
    
    int nref;
  };
      
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
    
    void addCell(DCell* c,
		 std::vector<DVertex*> vs, int parent_id);
    void removeCell(DCell* c);

    void imp(Mesh& mesh);
    void exp(Mesh& mesh);
    void number();

    void bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1);

    void bisectMarked(std::vector<bool> marked_ids);

    DCell* opposite(DCell* dcell, DVertex* v1, DVertex* v2);

    void propagate_naive(std::vector<uint>& propagated, bool& empty);

    void propagate_hypercube(std::vector<Propagation>& propagated, bool& empty);
    
    std::set<DVertex *> vertices;
    std::list<DCell *> cells;

    CellType* cell_type;
    uint d;

    // Start offset for new global id
    uint _start_offset, _max;

    // Propagation buffer
    std::vector<Propagation> propagate, type2, type3;

    // Map between global number of boundary vertex to vertex
    _map<uint, DVertex*> bc_dvs;
    std::map<EdgeKey, DVertex*> ref_edge;

    // Comparison operator for index/value pairs
    struct less_pair : public std::binary_function<std::pair<uint, prop_edge>,
						   std::pair<uint, prop_edge>, bool>
    {
      bool operator()(std::pair<uint, prop_edge> x, std::pair<uint, prop_edge> y)
      {
	return x.first < y.first;
      }
    };


    // Construct a edge id from given vertices
    inline EdgeKey edge_key(int id1, int id2) {
      dolfin_assert( id2 != id1);
      if(id2 == id1)
	error("Kaos id2 == id1");
      if(id2 < id1){
	EdgeKey key(id2,id1);    
	return key;
      }
      else {
	EdgeKey key(id1,id2);    
	return key;
      }      
    };

  };


}

#endif
