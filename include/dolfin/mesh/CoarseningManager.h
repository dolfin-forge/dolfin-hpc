// Copyright (C) 2013 Balthasar Reuter.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-03-25
// Last changed: 2013-04-03

#ifndef __COARSENING_MANAGER_H
#define __COARSENING_MANAGER_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/common/List.h>
#include <dolfin/mesh/MeshFunction.h>

namespace dolfin
{
  class Mesh;
  class Vertex;
  class DMesh;
  class DVertex;
  class DCell;

  /// Assists LocalMeshCoarsening by providing relevant information about mesh
  /// entities. 
  class CoarseningManager
  {
  public:
    class IndexMap;

    /// Create an empty instance, initialize it by calling init()
    CoarseningManager();

    /// Create a new instance and initialize it.
    explicit CoarseningManager(Mesh& mesh, MeshFunction<bool>& cell_marker, 
                               bool coarsen_boundary = false);

    ~CoarseningManager();

    /// Initialize the coarsening manager: build independent set, list of cells
    /// that are marked for coarsening and boundary information
    void init(Mesh& mesh, MeshFunction<bool>& cell_marker, 
              bool coarsen_boundary = false);

    /// Check if a vertex is forbidden, i. e. part of the independent set
    inline bool isForbiddenVertex(uint index)
    { return _forbidden_vertices.at(index); }

    /// Gives access to the dynamic mesh
    inline DMesh * dmesh()
    { return _dmesh; }

    /// Gives access to the dynamic mesh
    inline DMesh const * dmesh() const
    { return _dmesh; }

    /// Gives access to the list of cells for coarsening
    inline List< std::pair<DCell *, uint> >& cells_to_coarsen()
    { return _cells_to_coarsen; }

    /// Gives access to the list of cells for coarsening
    inline List< std::pair<DCell *, uint> > const & cells_to_coarsen() const
    { return _cells_to_coarsen; }

    /// Gives access to the list of cells that have been attempted to be 
    /// coarsened but failed because of missing data from other processes
    inline List<uint>& cells_to_request()
    { return _cells_to_request; }

    /// Gives access to the list of cells that have been attempted to be 
    /// coarsened but failed because of missing data from other processes
    inline List<uint> const & cells_to_request() const
    { return _cells_to_request; }

    /// Gives access to the list of vertices that have been attempted to be 
    /// deleted but failed because of missing data from other processes
    inline List<uint>& vertices_to_request()
    { return _vertices_to_request; }

    /// Gives access to the list of vertices that have been attempted to be 
    /// deleted but failed because of missing data from other processes
    inline List<uint> const & vertices_to_request() const
    { return _vertices_to_request; }

    /// Migrates cells according to request list
    bool migrate(uint num_cells_coarsened);

  private:
    /// Dynamic mesh used instead of the regular mesh object during coarsening
    DMesh* _dmesh;

    /// Indicator for independent set of vertices
    Array<bool> _forbidden_vertices;

    /// List of cells to coarsen
    List< std::pair<DCell *, uint> > _cells_to_coarsen;

    /// Original number of cells in the mesh
    uint _orig_num_cells;

    /// Original number of vertices in the mesh
    uint _orig_num_vertices;

    /// List of cells that need neighboring cells from other processes
    List<uint> _cells_to_request;

    /// List of vertices that need neighboring cells from other processes
    List<uint> _vertices_to_request;

    /// Maximum number of cells that have been coarsened in the last iteration
    uint _global_max_coarsened_cells;

    /// Maximum number of cells that remain for coarsening in the last iteration
    uint _global_max_remaining_cells;

    /// Number of migrations performed
    uint _migrations;

    /// Number of load balances performed
    uint _load_balances;

    /// performs part of the initialization which is also used by migrate().
    ///
    /// cell_marker is templated, because migrate() needs it for uint while
    /// init() provides bool. This way a conversion of the MeshFunction can
    /// be avoided.
    template<typename T>
    void initCommon(Mesh& mesh, MeshFunction<T>& cell_marker,
                    MeshFunction<uint> * attempt_count );

    /// Extracts an independent set of vertices by a greedy algorithm and
    /// stores it in _forbidden_vertices
    void findIndependentSet(Mesh& mesh, bool coarsen_boundary);

    /// Checks if the given vertex has neighbours that are part of the
    /// independent set and returns false in that case, otherwise true.
    bool isIndependentVertex(Vertex& v);

    /// Find all cells that are marked for coarsening and put them into a list
    template<typename T>
    void findCellsToCoarsen(MeshFunction<T>& cell_marker,
                            MeshFunction<uint> * attempt_count);

    /// Remove erased cells from the coarsening list
    void removeErasedCellsFromCoarseningList();

    /// Build Arrays of MeshFunctions that are required for the distribution
    /// step
    void buildMFArrays(Mesh& mesh, Array<int>& old2new_cells,
                       Array<int>& old2new_vertices,
                       Array< std::pair< MeshFunction<uint>*, 
                       MeshFunction<uint>* > >& cell_functions,
                       Array< std::pair< MeshFunction<double>*,
                       MeshFunction<double>*> >& vertex_functions);

    /// Cleanup the MeshFunction-Arrays after the distribution
    void cleanupMFArrays(Array< std::pair< MeshFunction<uint>*, 
                         MeshFunction<uint>* > >& cell_functions,
                         Array< std::pair< MeshFunction<double>*,
                         MeshFunction<double>*> >& vertex_functions);

    /// Update the independent set with the received values
    void updateIndependentSet(Mesh& mesh, 
                              MeshFunction<double>& forbidden_vertices_new);

    /// exchange requests for vertices
    void exchangeRequests(Mesh& mesh, Array<int>& old2new_cells,
                          Array<int>& old2new_vertices,
                          uint max_num_requested_vertices,
                          MeshFunction<uint> *& partitions);
  };
}

#endif 
