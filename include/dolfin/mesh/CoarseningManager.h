// Copyright (C) 2013 Balthasar Reuter.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-03-25
// Last changed: 2013-04-03

#ifndef __DOLFIN_COARSENING_MANAGER_H
#define __DOLFIN_COARSENING_MANAGER_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/common/List.h>
#include <dolfin/mesh/MeshFunction.h>

namespace dolfin
{

class Mesh;
class Vertex;
class DMesh;
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
  void init(Mesh& mesh, MeshFunction<bool>& cell_marker, bool coarsen_boundary =
                false);

  /// Check if a vertex is forbidden, i. e. part of the independent set
  inline bool isForbiddenVertex(uint index)
  {
    return forbidden_vertices_.at(index);
  }

  /// Gives access to the dynamic mesh
  inline DMesh * dmesh()
  {
    return dmesh_;
  }

  /// Gives access to the dynamic mesh
  inline DMesh const * dmesh() const
  {
    return dmesh_;
  }

  /// Gives access to the list of cells for coarsening
  inline List<std::pair<DCell *, uint> >& cells_to_coarsen()
  {
    return cells_to_coarsen_;
  }

  /// Gives access to the list of cells for coarsening
  inline List<std::pair<DCell *, uint> > const & cells_to_coarsen() const
  {
    return cells_to_coarsen_;
  }

  /// Gives access to the list of cells that have been attempted to be
  /// coarsened but failed because of missing data from other processes
  inline List<uint>& cells_to_request()
  {
    return cells_to_request_;
  }

  /// Gives access to the list of cells that have been attempted to be
  /// coarsened but failed because of missing data from other processes
  inline List<uint> const & cells_to_request() const
  {
    return cells_to_request_;
  }

  /// Gives access to the list of vertices that have been attempted to be
  /// deleted but failed because of missing data from other processes
  inline List<uint>& vertices_to_request()
  {
    return vertices_to_request_;
  }

  /// Gives access to the list of vertices that have been attempted to be
  /// deleted but failed because of missing data from other processes
  inline List<uint> const & vertices_to_request() const
  {
    return vertices_to_request_;
  }

  /// Migrates cells according to request list
  bool migrate(uint num_cells_coarsened);

  /// Returns the acceptable mesh quality threshold
  real qualityThreshold() const
  {
    return quality_threshold_;
  }

private:

  /// Dynamic mesh used instead of the regular mesh object during coarsening
  DMesh* dmesh_;

  /// Indicator for independent set of vertices
  Array<bool> forbidden_vertices_;

  /// List of cells to coarsen
  List<std::pair<DCell *, uint> > cells_to_coarsen_;

  /// Original number of cells in the mesh
  uint orig_num_cells_;

  /// Original number of vertices in the mesh
  uint orig_num_vertices_;

  /// List of cells that need neighboring cells from other processes
  List<uint> cells_to_request_;

  /// List of vertices that need neighboring cells from other processes
  List<uint> vertices_to_request_;

  /// Number of cells that have been migrated away by this process in the
  /// previous migration
  uint num_migrated_cells_;

  /// Quality threshold for the checkMesh routine
  real quality_threshold_;

  /// performs part of the initialization which is also used by migrate().
  ///
  /// attempt count gives the number of precious attempts, with the number
  /// being one larger than the actual number of attempts. This way
  /// 0 corresponds to not being marked for coarsening, 1 being marked for
  /// coarsening without prior attempts and >1 marked for coarsening with
  /// the value being the number of prior attempts minus one.
  void initCommon(Mesh& mesh, MeshFunction<uint> * attempt_count);

  /// Extracts an independent set of vertices by a greedy algorithm and
  /// stores it in _forbidden_vertices
  void findIndependentSet(Mesh& mesh, bool coarsen_boundary);

  /// Checks if the given vertex has neighbours that are part of the
  /// independent set and returns false in that case, otherwise true.
  bool isIndependentVertex(Vertex& v);

  /// Find all cells that are marked for coarsening and put them into a list
  /// together with the number of previous attempts
  void findCellsToCoarsen(MeshFunction<uint> * attempt_count);

  /// Remove erased cells from the coarsening list
  void removeErasedCellsFromCoarseningList();

  /// Build Arrays of MeshFunctions that are required for the distribution
  /// step
  void buildMFArrays(
      Mesh& mesh,
      Array<int>& old2new_cells,
      Array<int>& old2new_vertices,
      Array<std::pair<MeshFunction<uint>*, MeshFunction<uint>*> >& cell_functions,
      Array<std::pair<MeshFunction<real>*, MeshFunction<real>*> >& vertex_functions);

  /// Cleanup the MeshFunction-Arrays after the distribution
  void cleanupMFArrays(
      Array<std::pair<MeshFunction<uint>*, MeshFunction<uint>*> >& cell_functions,
      Array<std::pair<MeshFunction<real>*, MeshFunction<real>*> >& vertex_functions);

  /// Update the independent set with the received values
  void updateIndependentSet(Mesh& mesh,
                            MeshFunction<real>& forbidden_vertices_new);

  /// exchange requests for vertices
  void exchangeRequests(Mesh& mesh, Array<int>& old2new_cells,
                        Array<int>& old2new_vertices,
                        uint max_num_requested_vertices,
                        MeshFunction<uint> *& partitions);
};
}

#endif 
