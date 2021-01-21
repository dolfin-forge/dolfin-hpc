// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFC_CELL_H
#define __DOLFIN_UFC_CELL_H

#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/ufc/ufc.h>

namespace dolfin
{

/// This class is simple wrapper for a UFC cell and provides
/// a layer between a DOLFIN cell and a UFC cell.

class UFCCell : public ufc::cell
{
public:
  /// Create empty UFC cell
  UFCCell();

  /// Create UFC cell from DOLFIN cell
  UFCCell( Cell & dolfin_cell );

  /// Copy constructor
  UFCCell( UFCCell const & other );

  /// Destructor
  ~UFCCell();

public:
  /// Dereference operator, returns a reference to the underlying Cell
  inline auto operator*() const -> Cell const &;

  // Initialize UFC cell data
  void init( Cell & cell );

  // Update cell entities to global indices and coordinates
  void update( Cell & cell );

  ///
  static auto shape( CellType::Type type ) -> ufc::shape;

  // Number of cell vertices
  uint num_vertices { 0 };

  // coordinates
  std::vector< double > coordinates;

private:
  // access to dolfin cell
  Cell const * cell_ { nullptr };

  // Clear UFC cell data
  void clear();
};

//--- INLINES -----------------------------------------------------------------

inline UFCCell::UFCCell()
  : ufc::cell()
{
}

//-----------------------------------------------------------------------------

inline UFCCell::UFCCell( Cell & dolfin_cell )
  : ufc::cell()
  , cell_( &dolfin_cell )
{
  init( dolfin_cell );
}

//-----------------------------------------------------------------------------

inline UFCCell::UFCCell( UFCCell const & other )
  : ufc::cell()
  , cell_( other.cell_ )
{
  if ( cell_ != nullptr )
    init( *const_cast< Cell * >( cell_ ) );
}

//-----------------------------------------------------------------------------

inline UFCCell::~UFCCell()
{
  clear();
}

//-----------------------------------------------------------------------------

inline auto UFCCell::operator*() const -> Cell const &
{
  return *cell_;
}

//-----------------------------------------------------------------------------

inline auto UFCCell::shape( CellType::Type type ) -> ufc::shape
{
  switch ( type )
  {
    case CellType::interval:
      return ufc::shape::interval;
      break;
    case CellType::triangle:
      return ufc::shape::triangle;
      break;
    case CellType::tetrahedron:
      return ufc::shape::tetrahedron;
      break;
    case CellType::quadrilateral:
      return ufc::shape::quadrilateral;
      break;
    case CellType::hexahedron:
      return ufc::shape::hexahedron;
      break;
    default:
      error( "UFCCell : unknown cell type." );
      break;
  }
  return ufc::shape::interval;
}

//-----------------------------------------------------------------------------

inline void UFCCell::init( Cell & cell )
{
  // Clear old data
  clear();

  // Update dolfin cell pointer
  this->cell_ = &cell;

  // Set cell shape
  cell_shape = shape( cell.type() );

  // Set topological dimension
  topological_dimension = cell.dim();

  // Set geometric dimension
  geometric_dimension = cell.mesh().geometry_dimension();

  // Set entity indices
  entity_indices.resize( topological_dimension + 1 );

  // In any case store topological data in object
  for ( uint d = 0; d < topological_dimension; ++d )
  {
    entity_indices[d] = cell.entities( d );
  }

  // Cell index (short-cut for entity_indices[topological_dimension][0])
  entity_indices[topological_dimension] = { cell.global_index() };
  index = entity_indices[topological_dimension][0];

  // FIXME the following three are only set with placeholders
  local_facet = -1;
  orientation = -1;
  mesh_identifier = -1;

  //
  num_vertices = cell.num_entities( 0 );

  /// Set vertex coordinates
  Array< uint > const & vertices = cell.entities( 0 );
  coordinates.resize( num_vertices * geometric_dimension );

  for ( uint i = 0; i < num_vertices; ++i )
  {
    double const * coords = cell.mesh().geometry().x( vertices[i] );

    for ( uint c = 0; c < geometric_dimension; ++c )
      coordinates[i * geometric_dimension + c] = coords[c];
  }
}

//-----------------------------------------------------------------------------

inline void UFCCell::clear()
{
  entity_indices.clear();
  coordinates.clear();

  cell_shape            = ufc::shape::interval;
  topological_dimension = 0;
  geometric_dimension   = 0;
}

//-----------------------------------------------------------------------------

inline void UFCCell::update( Cell & cell )
{
  // Update dolfin cell pointer
  this->cell_ = &cell;

#if ENABLE_P1_OPTIMIZATIONS
  cell.get_global_entities( 0, entity_indices[0].data() );
#else
  cell.get_global_entities( entity_indices );
#endif
  entity_indices[topological_dimension][0] = cell.global_index();

  // Cell index (short-cut for entity_indices[topological_dimension][0])
  index = entity_indices[topological_dimension][0];

  // /// Set vertex coordinates
  Array< uint > const & vertices = cell.entities( 0 );
  for ( uint i = 0; i < num_vertices; ++i )
  {
    double const * coords = cell.mesh().geometry().x( vertices[i] );

    for ( uint c = 0; c < geometric_dimension; ++c )
      coordinates[i * geometric_dimension + c] = coords[c];
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_CELL_H */
