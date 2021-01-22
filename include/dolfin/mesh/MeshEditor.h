// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_EDITOR_H
#define __DOLFIN_MESH_EDITOR_H

#include <dolfin/common/types.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/CellType.h>

namespace dolfin
{

class Mesh;
class Connectivity;

//-----------------------------------------------------------------------------

/// A simple mesh editor for creating meshes.

class MeshEditor
{

public:
  /// Constructor for serial meshes given a cell type and default space
  MeshEditor( Mesh & mesh, CellType const & ctype );

  /// Constructor for serial meshes given a cell type and default space
  MeshEditor( Mesh & mesh, CellType const & ctype, Comm & comm );

  /// Constructor for serial meshes given a cell type and a space definition
  MeshEditor( Mesh & mesh, CellType const & ctype, Space const & space );

  /// Constructor for meshes given a cell type and a space definition
  MeshEditor( Mesh &           mesh,
              CellType const & ctype,
              Space const &    space,
              Comm &           comm );

  /// Constructor for serial meshes with type of cell from factory function
  MeshEditor( Mesh & mesh, CellType::Type cell_type, size_t gdim );

  /// Constructor for serial meshes with type of cell from factory function
  MeshEditor( Mesh & mesh, CellType::Type cell_type, size_t gdim, Comm & comm );

  /// Constructor using already initialized mesh.
  MeshEditor( Mesh & mesh );

  /// Destructor
  ~MeshEditor();

  ///
  auto mesh() const -> Mesh &
  {
    return mesh_;
  }

  //--- VERTICES --------------------------------------------------------------

  /// Specify number of vertices
  /// Optionally specify the global number of vertices for a distributed mesh.
  /// If the topology is not distributed, any value different than zero or the
  /// number of local vertices will trigger an error.
  auto init_vertices( size_t num_local, size_t num_global = 0 ) -> void;

  /// Add vertex v at given coordinates x
  auto add_vertex( size_t v, real const * x ) -> void;

  /// Return current vertex count
  auto current_vertex() const -> size_t;

  //--- CELLS -----------------------------------------------------------------

  /// Specify number of cells
  /// Optionally specify the global number of cells for a distributed mesh.
  /// If the topology is not distributed, any value different than zero or the
  /// number of local cells will trigger an error.
  auto init_cells( size_t num_local, size_t num_global = 0 ) -> void;

  /// Add cell with given vertices
  auto add_cell( size_t c, size_t const * v ) -> void;

  /// Return current cell count
  auto current_cell() const -> size_t;

  //---------------------------------------------------------------------------

  /// Close mesh, finish editing
  auto close() -> void;

  //---------------------------------------------------------------------------

private:
  /// Open mesh of given cell type and geometrical dimension
  auto init( Mesh &           mesh,
             CellType const & ctype,
             Space const &    space,
             Comm &           comm ) -> void;

  // Clear all data
  auto clear() -> void;

  // Mesh
  Mesh & mesh_;

  // Cell connectivity to vertices
  Connectivity * cell_vertices_;

  // Topological dimension
  size_t tdim_;

  // Geometrical (Euclidean) dimension
  size_t gdim_;

  // Number of vertices
  size_t num_vertices_;

  // Number of cells
  size_t num_cells_;

  // Next available vertex
  size_t vertex_index_;

  // Next available cell
  size_t cell_index_;

  //
  bool open_;
};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_EDITOR_H */
