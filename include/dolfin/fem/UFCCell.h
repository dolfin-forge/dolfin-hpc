// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFC_CELL_H
#define __DOLFIN_UFC_CELL_H

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshDistributedData.h>

#include <ufc.h>

namespace dolfin
{

/// This class is simple wrapper for a UFC cell and provides
/// a layer between a DOLFIN cell and a UFC cell.

class UFCCell : public ufc::cell
{
  friend class UFCCellIterator;

public:

  /// Create empty UFC cell
  UFCCell() :
      ufc::cell(),
      cell_(NULL),
      num_vertices(0)
  {
  }

  /// Create UFC cell from DOLFIN cell
  UFCCell(Cell& dolfin_cell) :
      ufc::cell(),
      cell_(&dolfin_cell),
      num_vertices(0)
  {
    init(dolfin_cell);
  }

  /// Copy constructor
  UFCCell(UFCCell const& other) :
      ufc::cell(),
      cell_(other.cell_),
      num_vertices(0)
  {
    if (cell_ != NULL) init(*const_cast<Cell*>(cell_));
  }

  /// Destructor
  ~UFCCell()
  {
    clear();
  }

private:

  //
  Cell const * cell_;

public:

  /// Dereference operator, returns a reference to the underlying Cell
  inline Cell const& operator*() const { return *cell_; };

  // Number of cell vertices
  uint num_vertices;

  // Initialize UFC cell data
  void init(Cell& cell);

  // Update cell entities to global indices and coordinates
  void update(Cell& cell);

  ///
  static ufc::shape shape(CellType::Type type)
  {
    switch (type)
      {
      case CellType::interval:
        return ufc::interval;
        break;
      case CellType::triangle:
        return ufc::triangle;
        break;
      case CellType::tetrahedron:
        return ufc::tetrahedron;
        break;
      case CellType::quadrilateral:
        return ufc::quadrilateral;
        break;
      case CellType::hexahedron:
        return ufc::hexahedron;
        break;
      default:
        error("UFCCell : unknown cell type.");
        break;
      }
    return ufc::interval;
  }

private:

  // Clear UFC cell data
  void clear();

};

//--- INLINES -----------------------------------------------------------------

inline void UFCCell::init(Cell& cell)
{
  // Clear old data
  clear();

  // Update dolfin cell pointer
  this->cell_ = &cell;

  // Set cell shape
  cell_shape = shape(cell.type());

  //
  num_vertices = cell.num_entities(0);

  // Set topological dimension
  topological_dimension = cell.dim();

  // Set geometric dimension
  geometric_dimension = cell.mesh().geometry_dimension();

  // Set entity indices
  entity_indices = new uint*[topological_dimension + 1];
  entity_indices[topological_dimension] = new uint[1];
  entity_indices[topological_dimension][0] = cell.index();

  // Cell index (short-cut for entity_indices[topological_dimension][0])
  index = entity_indices[topological_dimension][0];

  // In any case store topological data in object
  for (uint d = 0; d < topological_dimension; ++d)
  {
    entity_indices[d] = new uint[cell.num_entities(d)];
    for (uint i = 0; i < cell.num_entities(d); ++i)
    {
      entity_indices[d][i] = (cell.entities(d))[i];
    }
  }

  /// Set vertex coordinates
  Array<uint> const & vertices = cell.entities(0);
  coordinates = new real*[num_vertices];
  for (uint i = 0; i < num_vertices; ++i)
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
  this->cell_ = &cell;

#if ENABLE_P1_OPTIMIZATIONS
  cell.get_global_entities(0, entity_indices[0]);
#else
  cell.get_global_entities(entity_indices);
#endif
  entity_indices[topological_dimension][0] = cell.index();

  // Cell index (short-cut for entity_indices[topological_dimension][0])
  index = entity_indices[topological_dimension][0];

  /// Set vertex coordinates
  Array<uint> const & vertices = cell.entities(0);
  for (uint i = 0; i < num_vertices; ++i)
  {
    coordinates[i] = cell.mesh().geometry().x(vertices[i]);
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_CELL_H */
