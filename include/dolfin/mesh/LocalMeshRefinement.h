// Copyright (C) 2006 Johan Hoffman.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-11-01

#ifndef __DOLFIN_LOCAL_MESH_REFINEMENT_H
#define __DOLFIN_LOCAL_MESH_REFINEMENT_H

#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

class Cell;
class Edge;
class Mesh;
class MeshEditor;

/**
 *  @class  LocalMeshRefinement
 *
 *  @brief  This class implements local mesh refinement for different mesh
 *          types.
 */

class LocalMeshRefinement
{
public:

  /// Refine simplicial mesh locally by edge bisection
  static void refineMeshByEdgeBisection(Mesh& mesh,
                                        MeshValues<bool, Cell>& cell_marker,
                                        bool refine_boundary = true,
                                        real tf = 0.0,
                                        real tb = 0.0,
                                        real ts = 0.0,
                                        bool balance = true);
private:

  /// Bisect edge of simplex cell
  static void bisectEdgeOfSimplexCell(Cell& cell, Edge& edge,
                                      uint& new_vertex,
                                      MeshEditor& editor,
                                      uint& current_cell);

};

}

#endif
