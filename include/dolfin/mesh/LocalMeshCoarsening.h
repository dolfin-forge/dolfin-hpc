// Copyright (C) 2006 Johan Hoffman.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Balthasar Reuter, 2013.
//
// First added:  2006-11-01
// Last changed: 2013-04-03

#ifndef __DOLFIN_LOCAL_MESH_COARSENING_H
#define __DOLFIN_LOCAL_MESH_COARSENING_H

#include <dolfin/mesh/MeshValues.h>

#include <utility>
#include <list>

namespace dolfin
{

  class Mesh;
  class Vertex;
  class Cell;
  class MeshEditor;
  class CoarseningManager;
  class DMesh;
  class DCell;
  class DVertex;

  /// This class implements local mesh coarsening for different mesh types.
  class LocalMeshCoarsening
  {
  public:
    /// Coarsen simplicial mesh locally by edge collapse 
    ///
    /// *Arguments*
    ///
    ///   cell_marker (MeshValues<bool, Cell>&)
    ///     Indicates cells for coarsening with true
    ///
    ///   coarsen_boundary (bool)
    ///     Enable or disable coarsening of boundary cells
    ///
    static void coarsenMeshByEdgeCollapse(MeshValues<bool, Cell>& cell_marker,
                                          bool coarsen_boundary = false); 

  private:
    /// Selects the shortest edge for coarsening in the specified cell, which 
    /// does not have two forbidden vertices.
    ///
    /// *Arguments*
    ///
    ///   c (DCell*)
    ///     The cell to be coarsened
    ///
    ///   manager (CoarseningManager&)
    ///     The Coarsening manager
    ///
    ///   vertices (DVertex **)
    ///     List of vertex indices, that gets filled with pointers to endpoints
    ///     of the found edge. Has to be able to hold at least two pointer.
    ///
    /// *Returns*
    ///
    ///   bool
    ///     True if a suitable edge has been found
    ///
    static bool selectEdge(DCell* c, CoarseningManager& manager, 
                           DVertex * vertices[]);

    /// Selects the vertex that will be deleted. If one of the vertices is
    /// forbidden the other is chosen else the vertices are chosen alternately
    /// in each attempt.
    ///
    /// *Arguments*
    ///
    ///   vertices (DVertex **)
    ///     Pointer to the two endpoints
    ///
    ///   manager (CoarseningManager&)
    ///     The Coarsening manager
    ///
    ///   attempts (uint)
    ///     Number of previous attempts to coarsen the cell
    ///   
    /// *Returns*
    ///
    ///   int
    ///     The index of the vertex chosen for deletion, -1 if no vertex chosen
    ///     because entities from other processes are needed
    ///
    static int selectVertex(DVertex * vertices[], CoarseningManager& manager,
                            uint attempts);

    /// Checks the cells adjacent to the removed vertex for wrong orientation
    /// and sufficient large ratio of volume to diameter (avoid stretched cells).
    ///
    /// *Arguments*
    ///
    ///   cells_to_regenerate (std::list<DCell *>&)
    ///     List of cells that has been changed
    ///
    ///   cells_to_regenerate_orient (std::vector<uint>&)
    ///     Previous orientations of the changed cells
    ///
    ///   quality_threshold (real)
    ///     Threshold for cell quality
    /// 
    /// *Returns*
    ///
    ///   bool
    ///     true if all cells are ok, false otherwise
    ///
    static bool checkMesh(std::list<DCell *>& cells_to_regenerate,
                          std::vector<uint>& cells_to_regenerate_orient,
                          real quality_threshold);

    /// Coarsen a selected cell by edge collapse. It is called from 
    /// coarsenMeshByEdgeCollapse().
    ///
    /// *Arguments*
    ///
    ///   manager (CoarseningManager&)
    ///     The CoarseningManager for meta data
    ///
    ///   cell_to_coarsen (DCell*)
    ///     Pointer to the cell that is chosen for coarsening
    ///
    ///   attempts (uint)
    ///     Number of previous attempts to coarsen this cell
    ///
    /// *Returns*
    ///
    ///   int
    ///
    ///     The number of cells deleted during the coarsening of the chosen cell
    ///     Special cases:
    ///
    ///     * 0 if the cell can't be coarsened (due to forbidden vertices etc.)
    ///
    ///     * -1 if coarsening failed due to the checkMesh failure
    ///
    ///     * -2 if coarsening failed due to missing entities from other processes
    ///
    static int coarsenCell(CoarseningManager& manager, DCell* cell_to_coarsen,
                           uint attempts);

  }; // end class LocalMeshCoarsening

}

#endif
