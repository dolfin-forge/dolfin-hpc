// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Hoffman 2007.
// Modified by Magnus Vikstrøm 2007.
// Modified by Garth N. Wells 2007.
// Modified by Balthasar Reuter, 2013.
//
// First added:  2006-05-08
// Last changed: 2013-03-22

#ifndef __MESH_H
#define __MESH_H

#include <string>
#include <dolfin/common/types.h>
#include <dolfin/common/Variable.h>
#include "ALEType.h"
#include "MeshTopology.h"
#include "MeshGeometry.h"
#include "CellType.h"
#include "MeshDistributedData.h"

namespace dolfin
{

  
  template <class T> class MeshFunction;
  class MeshData;

  /// A Mesh consists of a set of connected and numbered mesh entities.
  ///
  /// Both the representation and the interface are dimension-independent,
  /// but a concrete interface is also provided for standard named mesh
  /// entities:
  ///
  ///     Entity  Dimension  Codimension
  ///
  ///     Vertex      0           -
  ///     Edge        1           -
  ///     Face        2           -
  ///
  ///     Facet       -           1
  ///     Cell        -           0
  ///
  /// When working with mesh iterators, all entities and connectivity
  /// are precomputed automatically the first time an iterator is
  /// created over any given topological dimension or connectivity.
  ///
  /// Note that for efficiency, only entities of dimension zero
  /// (vertices) and entities of the maximal dimension (cells) exist
  /// when creating a Mesh. Other entities must be explicitly created
  /// by calling init(). For example, all edges in a mesh may be created
  /// by a call to mesh.init(1). Similarly, connectivities such as
  /// all edges connected to a given vertex must also be explicitly
  /// created (in this case by a call to mesh.init(0, 1)).

  
  class Mesh : public Variable
  {
  public:

    
    /// Create empty mesh
    Mesh();

    /// Copy constructor
    Mesh(const Mesh& mesh);

    /// Create mesh from data file
    Mesh(std::string filename);

    
    /// Destructor
    ~Mesh();

    /// Assignment
    const Mesh& operator=(const Mesh& mesh);

    /// Return number of vertices
    inline uint numVertices() const { return _topology.size(0); }

    /// Return number of edges
    inline uint numEdges() const { return _topology.size(1); }

    /// Return number of faces
    inline uint numFaces() const { return _topology.size(2); }

    /// Return number of facets
    inline uint numFacets() const { return _topology.size(_topology.dim() - 1); }

    /// Return number of cells
    inline uint numCells() const { return _topology.size(_topology.dim()); }

    /// Return coordinates of all vertices
    inline real* coordinates() { return _geometry.x(); }

    /// Return coordinates of all vertices
    inline const real* coordinates() const { return _geometry.x(); }

    /// Return connectivity for all cells
    inline uint* cells() { return _topology(_topology.dim(), 0)(); }

    /// Return connectivity for all cells
    inline const uint* cells() const { return _topology(_topology.dim(), 0)(); }

    /// Return number of entities of given topological dimension
    inline uint size(uint dim) const { return _topology.size(dim); }

    
    /// Return mesh topology (non-const version)
    inline MeshTopology& topology() { return _topology; }

    /// Return mesh topology (const version)
    inline const MeshTopology& topology() const { return _topology; }

    /// Return mesh geometry (non-const version)
    inline MeshGeometry& geometry() { return _geometry; }

    /// Return mesh geometry (const version)
    inline const MeshGeometry& geometry() const { return _geometry; }

    /// Return mesh distribution data
    inline MeshDistributedData& distdata() { return _distdata; }

    /// Return mesh distribution data (const version)
    const inline MeshDistributedData& distdata() const { return _distdata; }

    /// Return mesh data
    MeshData& data();

    /// Return mesh cell type
    inline CellType& type() { dolfin_assert(_cell_type); return *_cell_type; }

    /// Return mesh cell type
    inline const CellType& type() const { dolfin_assert(_cell_type); return *_cell_type; }

    /// Compute entities of given topological dimension and return number of entities
    uint init(uint dim);

    /// Compute connectivity between given pair of dimensions
    void init(uint d0, uint d1);

    /// Compute all entities and connectivity
    void init();

    /// Clear all mesh data
    void clear();

    /// Order all mesh entities (not needed if "mesh order entities" is set)
    void order();

    /// Return true iff topology is ordered according to the UFC numbering
    bool ordered() const;

    /// Refine mesh uniformly
    void refine();

    /// Refine mesh according to cells marked for refinement
    void refine(MeshFunction<bool>& cell_markers, bool refine_boundary = true,
		bool load_balance = true);

    /// Coarsen mesh uniformly
    void coarsen();

    /// Coarsen mesh according to cells marked for coarsening
    void coarsen(MeshFunction<bool>& cell_markers, bool coarsen_boundary = false);

    /// Move coordinates of mesh according to new boundary coordinates
    void move(Mesh& boundary, ALEType method=lagrange);

    
    /// Smooth mesh using Lagrangian mesh smoothing 
    void smooth();

    
    /// Partition mesh into num_processes partitions
    void partition(MeshFunction<uint>& partitions);

    /// Partition mesh into num_partitions partitions
    void partition(MeshFunction<uint>& partitions, uint num_partitions);

    /// Partition mesh into num_partitions = numProc with weights on vertices
    void partition(MeshFunction<uint>& partitions, MeshFunction<uint>& weight);

    /// Partition mesh into num_partitions = numProc
    void partition_geom(MeshFunction<uint>& partitions);

    
    /// Distribute a mesh according to a mesh function
    void distribute(MeshFunction<uint>& distribution);

    /// Distribute a mesh according to a mesh function and transfer marked cells
    void distribute(MeshFunction<uint>& distribution, 
                    MeshFunction<bool>& cell_markers,
                    MeshFunction<bool>& new_cell_markers);

    /// Distribute a mesh according to a mesh function and transfer cell functions
    ///
    /// cell_functions contains pairs as <old_function,new_function>
    void distribute(MeshFunction<uint>& distribution, 
                    Array< std::pair< MeshFunction<uint> *, MeshFunction<uint> * > 
                    >& cell_functions);

    /// Distribute a mesh according to a mesh function and transfer vertex functions
    ///
    /// vertex_functions contains pairs as <old_function,new_function>
    void distribute(MeshFunction<uint>& distribution, 
                    Array< std::pair< MeshFunction<double> *, MeshFunction<double> * > 
                    >& vertex_functions);

    /// Distribute a mesh according to a mesh function and transfer cell and
    /// vertex functions
    ///
    /// cell_functions contains pairs as <old_function,new_function>
    ///
    /// vertex_functions contains pairs as <old_function,new_function>
    void distribute(MeshFunction<uint>& distribution, 
                    Array< std::pair< MeshFunction<uint> *, MeshFunction<uint> * > 
                    >& cell_functions,
                    Array< std::pair< MeshFunction<double> *, MeshFunction<double> * > 
                    >& vertex_functions);

    /// Renumber mesh global numbering
    void renumber();

    /// Display mesh data
    void disp() const;

    
    /// Return a short desriptive string
    std::string str() const;

    /// Output
    friend LogStream& operator<< (LogStream& stream, const Mesh& mesh);

    
  private:

    // Friends
    friend class MeshEditor;
    friend class TopologyComputation;
    friend class MeshOrdering;
    friend class MPIMeshCommunicator;

    // Mesh topology
    MeshTopology _topology;

    // Mesh geometry
    MeshGeometry _geometry;

    // Auxiliary mesh data
    MeshData* _data;

    // Cell type
    CellType* _cell_type;

    /// Return true iff topology is ordered according to the UFC numbering
    bool _ordered;

    /// Distribued Mesh data
    MeshDistributedData _distdata;

  };

}

#endif
