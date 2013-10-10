// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2013.
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

#ifdef HAVE_LIBGEOM
namespace libgeom
{
  class Geometry;
}
#endif

namespace dolfin
{
  class Mesh;
  class DCell;
  class DVertex;
  class CellType;

  /// Dynamic mesh class for on-the-fly changes to the mesh. It is used by the
  /// recursive RivaraRefinement and the EdgeCollapse-MeshCoarsening.
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

    /// Remove a vertex
    ///
    /// Entity is just marked as deleted, but not yet erased
    void removeVertex(DVertex* v);

    /// Remove a cell
    ///
    /// Entity is just marked as deleted, but not yet erased
    void removeCell(DCell* c);

    /// Erase removed entities from datastructures
    ///
    /// removeVertex() and removeCell() only mark entitites for deletion but are
    /// not actually erased
    void eraseRemovedEntities();

    /// Find Vertex by its local id
    DVertex* getVertex(int local_id);

    /// Find Cell by its local id
    DCell* getCell(int local_id);

    /// Import an existing mesh
    ///
    /// The local numbering of cells and vertices from the mesh-object is
    /// adopted in the DMesh
    void imp(Mesh& mesh);

#ifdef HAVE_LIBGEOM
    /// Import an existing mesh along with geometric parameters
    void imp(Mesh& mesh, MeshFunction<int>& patch_id_list,
	     MeshFunction<float>& bnd_u,
	     MeshFunction<float>& bnd_v);
#endif
    
    /// Export to a regular mesh
    void exp(Mesh& mesh);
    
#ifdef HAVE_LIBGEOM
    /// Export to a regular mesh along with geometric parameters
    void exp(Mesh& mesh, MeshFunction<int>& patch_id_list,
	     MeshFunction<float>& bnd_u,
	     MeshFunction<float>& bnd_v);
#endif
    
    /// Export to a regular mesh but keep numbering in the DMesh
    ///
    /// An optional mapping between old and new indices is generated. The Arrays
    /// have to have at least the size of the old numbering (as it was at the 
    /// time of the import).
    void expKeepNumbering(Mesh& mesh, Array<int> * old2new_cells = 0, 
                          Array<int> * old2new_vertices = 0);
    
    /// Renumber mesh entities locally
    ///
    /// An optional mapping between old and new indices is generated. The Arrays
    /// have to have the size of the old numbering
    void number( Array<int> *old2new_cells=0, 
		 Array<int> *old2new_vertices=0 );
    
    /// Bisect cell dcell
    ///
    /// The edge for the bisection is given by hv0 and hv1 and hangv is the
    /// hanging node of the opposite cell
    void bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1);
    
#ifdef HAVE_LIBGEOM
    /// Bisect cell dcell with respect to geometry
    void bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1,
		libgeom::Geometry& geom);
#endif
    
    /// Bisect marked cells
    void bisectMarked(std::vector<bool> marked_ids);
    
#ifdef HAVE_LIBGEOM
    /// Bisect marked cells with respect to geometry
    void bisectMarked(std::vector<bool> marked_ids, libgeom::Geometry& geom);
#endif
    
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
