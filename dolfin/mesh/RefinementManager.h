// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-01-21
// Last changed: 2008-07-23

#ifndef __REFINEMENT_MANAGER_H
#define __REFINEMENT_MANAGER_H

#include <dolfin/common/Array.h>
#include <dolfin/common/types.h>
#include "MeshFunction.h"
#include <dolfin/main/MPI.h>
#include <map>
#include <set>

namespace dolfin
{
  class Mesh;
  class Cell;
  class Edge;
  class RefinementManager
  {
  public:

    RefinementManager();
    RefinementManager(Mesh& mesh);
    ~RefinementManager();

    void init(Mesh& mesh);

    // Map new vertices global number
    void map_new_vertices(Array<uint> shared_edge,
			  Mesh& oldmesh, Mesh& newmesh);

    // Add a new vertex inside the refinement manager
    void add_new_vertex(uint* edge,uint vertex, Mesh& mesh, bool shared);

    // 
    void mark_localboundary(Mesh& oldmesh,
			    MeshFunction<bool>& cell_marker,
			    uint& num_new_vertices,
			    uint& num_new_cells);

    // Propagate refinement to other processes
    void propagate_refinement(Array<uint> shared_edge,
			      Mesh& oldmesh, Mesh& newmesh,
			      Array<uint>& new_edge);


    // Add a new nonshared vertex
    inline void addVertex(uint vertex, Mesh& mesh)
    { if(MPI::numProcesses() > 1) add_new_vertex(0, vertex, mesh, false);}
    
    // Add a new vertex on a shared edge
    inline void addVertex(uint* edge, uint vertex, Mesh& mesh)
    { if(MPI::numProcesses() > 1) add_new_vertex(edge, vertex, mesh, true);}

    // Check if the edge lies between processors
    inline bool on_boundary(Edge& edge)
    { return (MPI::numProcesses() > 1 ? boundary_edge.get(edge) : false);}

    // Check if the cell has recvied a propagation
    inline bool forbidden_cell(Cell& cell)
    { return (MPI::numProcesses() > 1 ? cell_forbidden.get(cell) : false);}

    // Check if the edge are forbidden from propagated refinement
    inline bool forbidden_edge(Edge& edge)
    { return (MPI::numProcesses() > 1 ? edge_forbidden.get(edge) : false);}

    inline uint edge_refined(Cell& cell)
    { dolfin_assert(cell_refedge.count(cell.index())); 
      return cell_refedge[cell.index()]; }


  private:
    typedef std::pair<uint, uint> EdgeKey;
    
    // Set of boundary cell index
    std::set<uint> boundary_set;

    // Construct a key from edge vertices
    EdgeKey edge_key(uint id1, uint id2);

    std::map< EdgeKey, uint> new_global;
    std::map< EdgeKey, uint> new_vertex;
    
    std::map<EdgeKey, uint> edge_cell_map;
    std::map<EdgeKey, bool> refined_edge;

    std::map<uint, uint> cell_refedge;

    MeshFunction<bool> boundary_edge;

    MeshFunction<bool> cell_forbidden; 

    MeshFunction<bool> edge_forbidden;
    
    uint _start_offset;
    bool _refm_init;
  };
}

#endif
