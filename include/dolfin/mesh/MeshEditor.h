// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-16
// Last changed: 2006-11-17

#ifndef __MESH_EDITOR_H
#define __MESH_EDITOR_H

#include <dolfin/common/types.h>
#include "CellType.h"
#include <dolfin/common/Array.h>

namespace dolfin
{

class Mesh;
class Point;

/// A simple mesh editor for creating simplicial meshes in 1D, 2D and 3D.

class MeshEditor
{
public:
  
  /// Constructor for meshes with unique type of cell
  //TODO: Deprecate.
  MeshEditor(Mesh& mesh, CellType::Type type, uint tdim, uint gdim);

  /// Constructor for meshes with unique type of cell
  MeshEditor(Mesh& mesh, CellType::Type type, uint gdim);

  /// Destructor
  ~MeshEditor();

  /// Specify number of vertices
  void initVertices(uint num_vertices);

  /// Specify number of cells
  void initCells(uint num_cells);

  /// Add vertex v at given point p
  ///FIXME: Deprecate.
  void addVertex(uint v, Point const& p);

  /// Add vertex v at given coordinates x
  void addVertex(uint v, real const * x);

  /// Add vertex v at given coordinate x
  ///FIXME: Deprecate.
  void addVertex(uint v, real x);

  /// Add vertex v at given coordinate (x, y)
  ///FIXME: Deprecate.
  void addVertex(uint v, real x, real y);

  /// Add vertex v at given coordinate (x, y, z)
  ///FIXME: Deprecate.
  void addVertex(uint v, real x, real y, real z);

  /// Add cell with given vertices
  void addCell(uint c, Array<uint> const& v);

  /// Add cell (interval) with given vertices
  void addCell(uint c, uint v0, uint v1);

  /// Add cell (triangle) with given vertices
  void addCell(uint c, uint v0, uint v1, uint v2);

  /// Add cell (tetrahedron) with given vertices
  void addCell(uint c, uint v0, uint v1, uint v2, uint v3);

  /// Close mesh, finish editing
  void close();

private:

  /// Open mesh of given cell type and geometrical dimension
  void init(Mesh& mesh, CellType::Type type, uint gdim);

  // Add vertex, common part
  void addVertexCommon(uint v);

  // Add cell, common part
  void addCellCommon(uint v);

  // Clear all data
  void clear();

  // The mesh
  Mesh * const mesh_;

  // Topological dimension
  uint tdim_;

  // Geometrical (Euclidean) dimension
  uint gdim_;

  // Number of vertices
  uint num_vertices_;

  // Number of cells
  uint num_cells_;

  // Next available vertex
  uint vertex_index_;

  // Next available cell
  uint cell_index_;

  // Temporary storage for local cell data
  Array<uint> vertices;

};

}

#endif
