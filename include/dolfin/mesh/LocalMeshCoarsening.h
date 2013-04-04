// Copyright (C) 2006 Johan Hoffman.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Balthasar Reuter, 2013.
//
// First added:  2006-11-01
// Last changed: 2013-04-03

#ifndef __LOCAL_MESH_COARSENING_H
#define __LOCAL_MESH_COARSENING_H

#include <utility>

#include "MeshFunction.h"

namespace dolfin
{

  class Mesh;
  class Vertex;
  class Cell;
  class MeshEditor;
  class CoarseningManager;

  /// This class implements local mesh coarsening for different mesh types.
  class LocalMeshCoarsening
  {
  public:
    /// Coarsen simplicial mesh locally by edge collapse 
    ///
    /// *Arguments*
    ///
    ///   mesh (Mesh&)
    ///     The mesh to be coarsened
    ///
    ///   cell_marker (MeshFunction<bool>&)
    ///     Indicates cells for coarsening with true
    ///
    ///   coarsen_boundary (bool)
    ///     Enable or disable coarsening of boundary cells
    static void coarsenMeshByEdgeCollapse(Mesh& mesh, 
                                          MeshFunction<bool>& cell_marker,
                                          bool coarsen_boundary = false); 

  private:

    /// Selects an edge for coarsening in the specified cell, on which not both
    /// vertices are forbidden or on the boundary between two processes. If 
    /// coarsen_boundary is false, the edge may also not be on the global 
    /// domain boundary.
    ///
    /// *Arguments*
    ///
    ///   c (Cell&)
    ///     The cell which edges are tested
    ///
    ///   manager (CoarseningManager&)
    ///     The Coarsening manager
    ///
    /// *Returns*
    ///
    ///   int
    ///     The index of the edge or -1 if none is found (i.e. all vertices
    ///     are forbidden)
    static int selectEdge(Cell& c, CoarseningManager& manager);
    static bool selectEdge(Cell& c, CoarseningManager& manager, uint *vertices);

    /// Selects the vertex that will be deleted. If one of the vertices is
    /// forbidden the other is chosen. If one of the vertices is on a boundary
    /// (process or domain), the other is chosen. If both are allowed the 
    /// vertex with larger index is chosen.
    ///
    /// *Arguments*
    ///
    ///   e (Edge&)
    ///     The selected edge
    ///
    ///   manager (CoarseningManager&)
    ///     The Coarsening manager
    ///   
    ///   vertD (uint&)
    ///     The index (on the coarse mesh) of the vertex, that is selected for deletion.
    ///
    ///   vertR (uint&)
    ///     The index (on the coarse mesh) of the vertex, that is selected to remain, i. e.
    ///     on which the deleted vertex will be collapsed.
    static bool selectVertex(Edge& e, CoarseningManager& manager,
                             uint& vertD, uint& vertR);
    static int selectVertex(uint *vertices, CoarseningManager& manager);

    /// Regenerates the cells adjacent to the deleted vertex and inserts
    /// them in the MeshEditor.
    ///
    /// *Arguments*
    ///
    ///   mesh (Mesh const & mesh)
    ///     The fine mesh, from which the vertex_to_remove is removed.
    ///
    ///   editor (MeshEditor& editor)
    ///     The MeshEditor, that is used to build up the coarsened mesh and into
    ///     which the regenerated cells are inserted.
    ///
    ///   vertex_to_remove (Vertex&)
    ///     The vertex that has been selected for deletion and that will be replaced
    ///     by the vertex with index vertR
    ///
    ///   c_id (uint)
    ///     The cell_id in the coarse mesh that the first regenerated cell should be
    ///     assigned. The other regenerated cells should get following higher indices.
    ///
    ///   cells_to_remove (MeshFuncion<bool> const &)
    ///     Indicator that shows which cells in the fine mesh should be regenerated.
    ///
    ///   manager (CoarseningManager&)
    ///     The Coarsening manager
    static void regenerateCells(Mesh const & mesh, MeshEditor& editor, 
                                Vertex& vertex_to_remove, 
                                uint vertR, uint c_id, 
                                MeshFunction<bool> const & cells_to_remove,
                                CoarseningManager& manager);

    /// Checks the cells adjacent to the removed cell for wrong orientation
    /// and sufficient large ratio of volume to diameter (avoid stretched cells).
    ///
    /// *Arguments*
    ///
    ///   removed_cell (Cell&)
    ///     The removed cell, that has originally been chosen for coarsening
    ///
    ///   coarse_mesh (Mesh&)
    ///     The coarse mesh that will be checked
    ///
    ///   manager (CoarseningManager&)
    ///     The Coarsening manager
    /// 
    /// *Returns*
    ///   bool
    ///     true if all cells are ok, false otherwise
    static bool checkMesh(Vertex& removed_vertex, Mesh& coarse_mesh, 
                          CoarseningManager& manager);

    /// Coarsen a selected cell by edge collapse. Is called from 
    /// coarsenMeshByEdgeCollapse().
    ///
    /// *Arguments*
    ///
    ///   mesh (Mesh&)
    ///     The original mesh
    ///
    ///   coarse_mesh (Mesh&)
    ///     The coarsened mesh
    ///
    ///   manager (CoarseningManager&)
    ///     The CoarseningManager for meta data
    ///
    ///   cell_to_coarsen_id (uint)
    ///     Index of the cell to be coarsened
    ///
    /// *Returns*
    ///
    ///   std::pair<bool,bool>
    ///
    ///     The first value indicates wether the coarse mesh is ok, i.e. the 
    ///     return value of checkMesh(). If it's true coarsening was successful.
    ///     The second value indicates wether the mesh has been changed. Typically
    ///     the following return values will occur:
    ///
    ///     - (true,true): coarsening was successful, mesh quality is ok and the number 
    ///       of cells and vertices has been reduced during that process.
    ///     - (true,false): mesh quality is ok but the mesh has not been changed, usually
    ///       because all vertices of the selected cell are forbidden.
    ///     - (false,true): coarsening has been tried but the mesh quality is too bad
    ///       afterwards such that reverting the changes is recommended.
    ///
    static std::pair<bool,bool> coarsenCell(Mesh& mesh, Mesh& coarse_mesh, 
                                            CoarseningManager& manager,
                                            uint cell_to_coarsen_id);

  }; // end class LocalMeshCoarsening

}

#endif
