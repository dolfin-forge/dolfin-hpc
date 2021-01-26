// Copyright (C) 2013 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFC_REFERENCE_CELL_H
#define __DOLFIN_UFC_REFERENCE_CELL_H

#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/main/MPI.h>

#include <dolfin/ufc/ufc.h>

namespace dolfin
{

/// This class is a wrapper for a UFC cell and provides
/// a reference cell in DOLFIN.

class UFCReferenceCell: public ufc::cell
{
public:

  /// Create emtpy UFC cell
  UFCReferenceCell()
    : ufc::cell()
  {
  }

  /// Create UFC cell from DOLFIN cell
  UFCReferenceCell(Cell& cell) :
      ufc::cell(),
      num_vertices(0)
  {
    init(cell);
  }

  /// Destructor
  ~UFCReferenceCell() override
  {
    clear();
  }

  /// Initialize UFC cell data
  void init(Cell& cell)
  {
    // Clear old data
    clear();

    // Set cell shape
    switch (cell.type())
    {
      case CellType::interval:
        cell_shape = ufc::shape::interval;
        num_vertices = 2;
        break;
      case CellType::triangle:
        cell_shape = ufc::shape::triangle;
        num_vertices = 3;
        break;
      case CellType::tetrahedron:
        cell_shape = ufc::shape::tetrahedron;
        num_vertices = 4;
        break;
      default:
        error("Unknown cell type.");
        break;
    }

    // Set topological dimension
    topological_dimension = cell.dim();

    // Set geometric dimension
    geometric_dimension = cell.mesh().geometry_dimension();

    // Set entity indices
    entity_indices = new size_t*[topological_dimension + 1];
    entity_indices[topological_dimension] = new size_t[1];
    entity_indices[topological_dimension][0] = cell.index();

    // Two different cases
    if (MPI::size() == 1)
    {
      // Single process, pointer to mesh topological data
      for (size_t d = 0; d < topological_dimension; d++)
        entity_indices[d] = cell.entities(d).data();
    }
    else
    {
      // Parallel case, store topological data in object
      for (size_t d = 0; d < topological_dimension; d++)
      {
        entity_indices[d] = new size_t[cell.num_entities(d)];
        for (size_t i = 0; i < cell.num_entities(d); i++)
          entity_indices[d][i] = cell.entities(d)[i];
      }
    }

    /// Set vertex coordinates
    coordinates = new real*[num_vertices];
    for (size_t i = 0; i < num_vertices; i++)
      coordinates[i] = new real[3];

    switch (topological_dimension)
    {
      case 1:
        {
          coordinates[0][0] = 0.;
          coordinates[1][0] = 1.;
          break;
        }
      case 2:
        {
          coordinates[0][0] = 0.;
          coordinates[0][1] = 0.;
          coordinates[1][0] = 1.;
          coordinates[1][1] = 0.;
          coordinates[2][0] = 0.;
          coordinates[2][1] = 1.;
          break;
        }
      case 3:
        {
          coordinates[0][0] = 0.;
          coordinates[0][1] = 0.;
          coordinates[0][2] = 0.;
          coordinates[1][0] = 1.;
          coordinates[1][1] = 0.;
          coordinates[1][2] = 0.;
          coordinates[2][0] = 0.;
          coordinates[2][1] = 1.;
          coordinates[2][2] = 0.;
          coordinates[3][0] = 0.;
          coordinates[3][1] = 0.;
          coordinates[3][2] = 1.;
          break;
        }
    }
  }

  // Clear UFC cell data
  void clear()
  {
    if (entity_indices)
    {
      if (MPI::size() > 1)
      {
        for (size_t i = 0; i < (topological_dimension + 1); i++)
          delete[] entity_indices[i];
        delete[] entity_indices;
      }
      else
      {
        delete[] entity_indices[topological_dimension];
        delete[] entity_indices;
      }
    }
    entity_indices = 0;

    if (coordinates)
      delete[] coordinates;
    coordinates = 0;

    cell_shape = ufc::interval;
    topological_dimension = 0;
    geometric_dimension = 0;
  }

private:

  // Number of cell vertices
  size_t num_vertices{0};

};

}

#endif
