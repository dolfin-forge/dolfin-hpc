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

  /// Dynamic mesh class for on-the-fly changes to the mesh. It is used by the
  /// recursive RivaraRefinement.
  ///
  /// It provides import- and export-routines from and to the regular Mesh class
  /// of DOLFIN but cannot be used directly for MeshFunction or assembly
  /// routines.
  ///
  /// For mesh refinement it has methods for cell-based bisection.
  ///
  class DMesh
  {
  public:
    DMesh();
    ~DMesh();

    /// Edge Identifier containing ids of end vertices
    typedef std::pair<int, int> EdgeKey;    

    /// Edge data structure for propagation
    typedef struct __edge__ {
      uint mv;    //< global index of midpoint vertex
      uint v1;    //< global index of endpoint
      uint v2;    //< global index of endpoint
      uint owner; //< rank of owner
    } prop_edge;

    /// Pair datatype for propagation 
    typedef std::pair<uint, prop_edge> Propagation;

    /// Add a new vertex
    void addVertex(DVertex* v);
    
    /// Add a new cell with vertices vs and inside existing cell parent_id
    void addCell(DCell* c, std::vector<DVertex*> vs, int parent_id);

    /// Remove a cell
    void removeCell(DCell* c);

    /// Import an existing mesh
    void imp(Mesh& mesh);

    /// Export to a regular mesh
    void exp(Mesh& mesh);

    /// Renumber mesh entities locally
    void number();

    /// Bisect cell 
    ///
    /// TODO: what is hangv, hv0, hv1???
    void bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1);

    /// Bisect marked cells
    void bisectMarked(std::vector<bool> marked_ids);

    /// Get opposite cell with respect to vertices v1 and v2
    DCell* opposite(DCell* dcell, DVertex* v1, DVertex* v2);

    /// Propagate refinement
    ///
    /// TODO: what are the arguments???
    inline void propagate_refinement(std::vector<Propagation>& propagated,
                                     bool& empty)
    {
      if ( MPI::numProcesses() & ( MPI::numProcesses() - 1 ) )
        propagate_naive(propagated, empty); 
      else
        propagate_hypercube(propagated, empty); 
    }

    /// Naive refinement propagation with pairwise communication
    void propagate_naive(std::vector<Propagation>& propagated, bool& empty);

    /// Refinement propagation within hypercube
    void propagate_hypercube(std::vector<Propagation>& propagated, bool& empty);
    
    /// Vertices contained in the mesh
    std::set<DVertex *> vertices;

    /// Cells contained in the mesh
    std::list<DCell *> cells;

    /// CellType of mesh
    CellType* cell_type;

    /// Dimension of MeshTopology
    uint d;

    /// Maximum number of vertices in one process
    uint glb_max;

    /// enumeration salt for bisect
    uint _salt;

    /// Start offset for new global id
    uint _start_offset;

    /// Global maximum start offset
    uint _max;

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


    /// Construct an edge id from given vertices
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