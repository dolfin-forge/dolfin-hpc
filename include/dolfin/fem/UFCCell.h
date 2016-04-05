// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-03-01
// Last changed: 2007-03-13

#ifndef __DOLFIN_UFC_CELL_H
#define __DOLFIN_UFC_CELL_H

#include <ufc.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/main/MPI.h>

namespace dolfin
{

/// This class is simple wrapper for a UFC cell and provides
/// a layer between a DOLFIN cell and a UFC cell.

class UFCCell : public ufc::cell
{
public:

  /// Create emtpy UFC cell
  UFCCell() :
      ufc::cell(),
      cell(NULL),
      num_vertices_(0)
  {
  }

  /// Create UFC cell from DOLFIN cell
  UFCCell(Cell& dolfin_cell) :
      ufc::cell(),
      cell(&dolfin_cell),
      num_vertices_(0)
  {
    init(dolfin_cell);
  }

  /// Destructor
  ~UFCCell()
  {
    clear();
  }

  //
  Cell const * cell;

  // Initialize UFC cell data
  void init(Cell& cell);

  // Clear UFC cell data
  void clear();

  // Update cell entities to global indices and coordinates
  void update(Cell& cell);

private:

  // Number of cell vertices
  uint num_vertices_;

};

//--- INLINES -----------------------------------------------------------------

inline void UFCCell::init(Cell& cell)
{
  // Clear old data
  clear();

  // Update dolfin cell pointer
  this->cell = &cell;

  // Set cell shape
  switch (cell.type())
    {
    case CellType::interval:
      cell_shape = ufc::interval;
      break;
    case CellType::triangle:
      cell_shape = ufc::triangle;
      break;
    case CellType::tetrahedron:
      cell_shape = ufc::tetrahedron;
      break;
    default:
      error("Unknown cell type.");
      break;
    }
  num_vertices_ = cell.num_entities(0);

  // Set topological dimension
  topological_dimension = cell.dim();

  // Set geometric dimension
  geometric_dimension = cell.mesh().geometry().dim();

  // Set entity indices
  entity_indices = new uint*[topological_dimension + 1];
  entity_indices[topological_dimension] = new uint[1];
  entity_indices[topological_dimension][0] = cell.index();

  // Cell index (short-cut for entity_indices[topological_dimension][0])
  index = entity_indices[topological_dimension][0];

#if ENABLE_P1_OPTIMIZATIONS
  // Do no allocate edges/faces to make sure any invalid access triggers a
  // segmentation fault
  entity_indices[0] = new uint[cell.num_entities(0)];
  for (uint i = 0; i < cell.num_entities(0); ++i)
  {
    entity_indices[0][i] = (cell.entities(0))[i];
  }
  for (uint d = 1; d < topological_dimension; ++d)
  {
    entity_indices[d] = NULL;
  }
#else
  // In any case store topological data in object
  for (uint d = 0; d < topological_dimension; ++d)
  {
    entity_indices[d] = new uint[cell.num_entities(d)];
    for (uint i = 0; i < cell.num_entities(d); ++i)
    {
      entity_indices[d][i] = (cell.entities(d))[i];
    }
  }
#endif

  /// Set vertex coordinates
  uint* vertices = cell.entities(0);
  coordinates = new real*[num_vertices_];
  for (uint i = 0; i < num_vertices_; ++i)
  {
    coordinates[i] = cell.mesh().geometry().x(vertices[i]);
  }
}

//-----------------------------------------------------------------------------
inline void UFCCell::clear()
{
  if (entity_indices)
  {
    for (uint i = 0; i < (topological_dimension + 1); ++i)
    {
      delete[] entity_indices[i];
    }
    delete[] entity_indices;
  }
  entity_indices = 0;

  delete[] coordinates;
  coordinates = 0;

  cell_shape = ufc::interval;
  topological_dimension = 0;
  geometric_dimension = 0;
}

//-----------------------------------------------------------------------------
inline void UFCCell::update(Cell& cell)
{
  // Update dolfin cell pointer
  this->cell = &cell;

  // Cell index (short-cut for entity_indices[topological_dimension][0])
  index = cell.global_index();

#if ENABLE_P1_OPTIMIZATIONS
  cell.global_entities(0, entity_indices[0]);
#else
  for (uint d = 0; d < topological_dimension; ++d)
  {
    cell.global_entities(d, entity_indices[d]);
  }
#endif
  entity_indices[topological_dimension][0] = index;

  /// Set vertex coordinates
  uint const * vertices = cell.entities(0);
  for (uint i = 0; i < num_vertices_; ++i)
  {
    coordinates[i] = cell.mesh().geometry().x(vertices[i]);
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_CELL_H */
