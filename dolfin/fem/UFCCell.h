// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-03-01
// Last changed: 2007-03-13

#ifndef __UFC_CELL_H
#define __UFC_CELL_H

#include <ufc.h>

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
    UFCCell() : ufc::cell(), num_vertices(0) {}

    /// Create UFC cell from DOLFIN cell
    UFCCell(Cell& cell) : ufc::cell(), num_vertices(0)
    {
      init(cell);
    }

    /// Destructor
    ~UFCCell()
    {
      clear();
    }
    
    /// Initialize UFC cell data
    void init(Cell& cell)
    {
      // Clear old data
      clear();

      // Set cell shape
      switch ( cell.type() )
      {
      case CellType::interval:
        cell_shape = ufc::interval;
        num_vertices = 2;
        break;
      case CellType::triangle:
        cell_shape = ufc::triangle;
        num_vertices = 3;
        break;
      case CellType::tetrahedron:
        cell_shape = ufc::tetrahedron;
        num_vertices = 4;
        break;
      default:
        error("Unknown cell type.");
      }

      // Set topological dimension
      topological_dimension = cell.mesh().topology().dim();

      // Set geometric dimension
      geometric_dimension = cell.mesh().geometry().dim();
      
      // Set entity indices
      entity_indices = new uint*[topological_dimension + 1];
      entity_indices[topological_dimension] = new uint[1]; 
      entity_indices[topological_dimension][0] = cell.index();	
     
      // Two different cases
      if(MPI::numProcesses() == 1) {	
	// Single process, pointer to mesh topological data
	for (uint d = 0; d < topological_dimension; d++)
	  entity_indices[d] = cell.entities(d);
      }
      else {
	// Parallel case, store topological data in object
	for (uint d = 0; d < topological_dimension; d++) {
	  entity_indices[d] = new uint[cell.numEntities(d)];
	  for (uint i = 0; i < cell.numEntities(d); i++)
	    entity_indices[d][i] = (cell.entities(d))[i];
	}
      }

      /// Set vertex coordinates
      uint* vertices = cell.entities(0);
      coordinates = new real*[num_vertices];
      for (uint i = 0; i < num_vertices; i++)
        coordinates[i] = cell.mesh().geometry().x(vertices[i]);
    }

    // Clear UFC cell data
    void clear()
    {
      if ( entity_indices )
      {
	if(MPI::numProcesses() > 1) {
	  for (uint i = 0; i < (topological_dimension + 1); i++)
	    delete [] entity_indices[i];
	  delete[] entity_indices;
	}
	else {
	  delete [] entity_indices[topological_dimension];
	  delete [] entity_indices;
	}
      }
      entity_indices = 0;

      if ( coordinates )
        delete [] coordinates;
      coordinates = 0;

      cell_shape = ufc::interval;
      topological_dimension = 0;
      geometric_dimension = 0;
    }

    // Update cell entities and coordinates
    inline void update(Cell& cell, MeshDistributedData& distdata)
    {

      // Set entity indices
      if( MPI::numProcesses() == 1 ) {
	for (uint d = 0; d < topological_dimension; d++)
	  entity_indices[d] = cell.entities(d);
	entity_indices[topological_dimension][0] = cell.index();
      }
      else {
#if ENABLE_P1_OPTIMIZATIONS
	for(uint i = 0; i < cell.numEntities(0); i++)
	  entity_indices[0][i] = distdata.get_global( (cell.entities(0))[i], 0);
#else
	for (uint d = 0; d < topological_dimension; d++)
	  for(uint i = 0; i < cell.numEntities(d); i++)
	    entity_indices[d][i] = distdata.get_global( (cell.entities(d))[i], d);
#endif
	entity_indices[topological_dimension][0] = distdata.get_cell_global(cell.index());
      }
      
      /// Set vertex coordinates
      const uint* vertices = cell.entities(0);
      for (uint i = 0; i < num_vertices; i++)
        coordinates[i] = cell.mesh().geometry().x(vertices[i]);
    }

  private:
    
    // Number of cell vertices
    uint num_vertices;

  };

}

#endif
