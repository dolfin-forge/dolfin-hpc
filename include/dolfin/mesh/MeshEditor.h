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
#include <dolfin/mesh/Point.h>

namespace dolfin
{

class Mesh;
class MeshConnectivity;

/// A simple mesh editor for creating meshes in 1D, 2D and 3D.

class MeshEditor
{

public:

  /// Constructor for meshes given a cell type and a space definition
  MeshEditor(Mesh& mesh, CellType const& type, Space const& space);

  /// Constructor for meshes with unique type of cell from factory function
  MeshEditor(Mesh& mesh, CellType::Type cell_type, uint gdim);

  /// Constructor using already initialized mesh.
  MeshEditor(Mesh& mesh);

  /// Destructor
  ~MeshEditor();

  ///
  Mesh& mesh() const { return mesh_; }

  //--- VERTICES --------------------------------------------------------------

  /// Specify number of vertices
  /// Optionally specify the global number of vertices for a distributed mesh.
  /// If the topology is not distributed, any value different than zero or the
  /// number of local vertices will trigger an error.
  void init_vertices(uint num_local, uint num_global = 0);

  /// Add vertex v at given coordinates x
  void add_vertex(uint v, real const * x);

  /// Return current vertex count
  uint current_vertex() const;

  //--- CELLS -----------------------------------------------------------------

  /// Specify number of cells
  /// Optionally specify the global number of cells for a distributed mesh.
  /// If the topology is not distributed, any value different than zero or the
  /// number of local cells will trigger an error.
  void init_cells(uint num_local, uint num_global = 0);

  /// Add cell with given vertices
  void add_cell(uint c, uint const * v);

  /// Return current cell count
  uint current_cell() const;

  //---------------------------------------------------------------------------

  /// Close mesh, finish editing
  void close();

private:

  /// Open mesh of given cell type and geometrical dimension
  void init(Mesh& mesh, CellType const& type, Space const& space);

  // Clear all data
  void clear();

  // Mesh
  Mesh& mesh_;

  // Cell connectivity to vertices
  MeshConnectivity * cell_vertices_;

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
