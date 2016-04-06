// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-16
// Last changed: 2006-11-17

#ifndef __DOLFIN_MESH_EDITOR_H
#define __DOLFIN_MESH_EDITOR_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/mesh/CellType.h>

namespace dolfin
{

class Mesh;
class Point;

/// A simple mesh editor for creating simplicial meshes in 1D, 2D and 3D.

class MeshEditor
{
public:

  /// Constructor for meshes with unique type of cell
  MeshEditor(Mesh& mesh, CellType const& cell_type, uint gdim);

  /// Constructor for meshes with unique type of cell from factory function
  MeshEditor(Mesh& mesh, CellType::Type type, uint gdim);

  /// Destructor
  ~MeshEditor();

  /// Specify number of vertices
  void init_vertices(uint num_vertices);

  /// Specify number of cells
  void init_cells(uint num_cells);

  /// Add vertex v at given coordinates x
  void add_vertex(uint v, real const * x);

  /// Add cell with given vertices
  void add_cell(uint c, uint const * v);

  /// Close mesh, finish editing
  void close();

  /// Return current vertex count
  uint current_vertex() const;

  /// Return current cell count
  uint current_cell() const;

private:

  /// Open mesh of given cell type and geometrical dimension
  void init(Mesh& mesh, CellType const& type, uint gdim);

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

  //
  bool open_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_EDITOR_H */
