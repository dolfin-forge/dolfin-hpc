// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFC_CELL_H
#define __DOLFIN_UFC_CELL_H

#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/MeshDistributedData.h>

#include <utility>

#include <ufc.h>
#include <utility>

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
  UFCCell( Cell const & dolfin_cell );

  /// Copy constructor
  UFCCell( UFCCell const & copy );

  /// Move constructor
  UFCCell( UFCCell && move );

  /// Destructor
  ~UFCCell();

  auto operator=( UFCCell const & copy ) -> UFCCell &;
  auto operator=( UFCCell && move ) -> UFCCell &;


public:
  /// Dereference operator, returns a reference to the underlying Cell
  inline auto operator*() const -> Cell const &;

  // Initialize UFC cell data
  auto init( Cell const & cell ) -> void;

  // Update cell entities to global indices and coordinates
  auto update( Cell const & cell ) -> void;

  ///
  static auto shape( CellType::Type type ) -> ufc::shape;

  // Number of cell vertices
  size_t num_vertices { 0 };

  // coordinates
  std::vector< double > coordinates;

private:
  // access to dolfin cell
  Cell const * cell_ { nullptr };

  // Clear UFC cell data
  auto clear() -> void;
};

//--- INLINES -----------------------------------------------------------------

inline UFCCell::UFCCell()
  : ufc::cell()
{
  cell_shape            = ufc::shape::interval;
  topological_dimension = DOLFIN_SIZE_T_UNDEF;
  geometric_dimension   = DOLFIN_SIZE_T_UNDEF;
  entity_indices        = {};
  index                 = DOLFIN_SIZE_T_UNDEF;
  local_facet           = DOLFIN_INT_UNDEF;
  orientation           = DOLFIN_INT_UNDEF;
  mesh_identifier       = DOLFIN_INT_UNDEF;
}

//-----------------------------------------------------------------------------

inline UFCCell::UFCCell( Cell const & dolfin_cell )
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
  {
    init( *cell_ );
  }
  else
  {
    cell_shape            = ufc::shape::interval;
    topological_dimension = DOLFIN_SIZE_T_UNDEF;
    geometric_dimension   = DOLFIN_SIZE_T_UNDEF;
    entity_indices        = {};
    index                 = DOLFIN_SIZE_T_UNDEF;
    local_facet           = DOLFIN_INT_UNDEF;
    orientation           = DOLFIN_INT_UNDEF;
    mesh_identifier       = DOLFIN_INT_UNDEF;
  }
}

//-----------------------------------------------------------------------------

inline UFCCell::UFCCell( UFCCell && move )
  : ufc::cell( std::move( move ) )
  , num_vertices( std::exchange( move.num_vertices, 0 ) )
  , coordinates( std::move( move.coordinates ) )
  , cell_( std::exchange( move.cell_, nullptr ) )
{
}

//-----------------------------------------------------------------------------

inline UFCCell::~UFCCell()
{
  entity_indices.clear();
  coordinates.clear();

  cell_shape            = ufc::shape::interval;
  topological_dimension = 0;
  geometric_dimension   = 0;
}

//-----------------------------------------------------------------------------

inline auto UFCCell::operator=( UFCCell const & copy ) -> UFCCell &
{
  this->cell_shape            = copy.cell_shape;
  this->topological_dimension = copy.topological_dimension;
  this->geometric_dimension   = copy.geometric_dimension;

  this->entity_indices.resize( copy.entity_indices.size() );
  for ( size_t i = 0; i < this->entity_indices.size(); ++i )
  {
    this->entity_indices[i] = copy.entity_indices[i];
  }

  this->index           = copy.index;
  this->local_facet     = copy.local_facet;
  this->orientation     = copy.orientation;
  this->mesh_identifier = copy.mesh_identifier;

  this->num_vertices = copy.num_vertices;
  this->coordinates  = copy.coordinates;
  this->cell_        = copy.cell_;

  return *this;
}

//-----------------------------------------------------------------------------

inline auto UFCCell::operator=( UFCCell && move ) -> UFCCell &
{
  this->cell_shape = std::exchange( move.cell_shape, ufc::shape::interval );
  this->topological_dimension = std::exchange( move.topological_dimension, 0 );
  this->geometric_dimension   = std::exchange( move.geometric_dimension, 0 );
  this->entity_indices        = std::move( move.entity_indices );
  this->index                 = std::exchange( move.index, 0 );
  this->local_facet           = std::exchange( move.local_facet, 0 );
  this->orientation           = std::exchange( move.orientation, 0 );
  this->mesh_identifier       = std::exchange( move.mesh_identifier, 0 );

  this->num_vertices = std::exchange( move.num_vertices, 0 );
  this->coordinates  = std::move( move.coordinates );
  this->cell_        = std::exchange( move.cell_, nullptr );

  return *this;
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

inline auto UFCCell::init( Cell const & cell ) -> void
{
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
  for ( size_t d = 0; d < topological_dimension; ++d )
  {
    entity_indices[d] = cell.entities( d );
  }

  // Cell index (short-cut for entity_indices[topological_dimension][0])
  entity_indices[topological_dimension] = { cell.global_index() };
  index = entity_indices[topological_dimension][0];

  // FIXME the following three are only set with placeholders
  local_facet = -1;
  orientation = 0;
  mesh_identifier = -1;

  //
  num_vertices = cell.num_entities( 0 );

  /// Set vertex coordinates
  std::vector< size_t > const & vertices = cell.entities( 0 );
  coordinates.resize( num_vertices * geometric_dimension );

  for ( size_t i = 0; i < num_vertices; ++i )
  {
    double const * coords = cell.mesh().geometry().x( vertices[i] );

    for ( size_t c = 0; c < geometric_dimension; ++c )
      coordinates[i * geometric_dimension + c] = coords[c];
  }
}

//-----------------------------------------------------------------------------

inline auto UFCCell::update( Cell const & cell ) -> void
{
  // Update dolfin cell pointer
  this->cell_ = &cell;

#if ENABLE_P1_OPTIMIZATIONS
  cell.get_global_entities( 0, entity_indices[0] );
#else
  cell.get_global_entities( entity_indices );
#endif
  entity_indices[topological_dimension][0] = cell.global_index();

  // Cell index (short-cut for entity_indices[topological_dimension][0])
  index = entity_indices[topological_dimension][0];

  // /// Set vertex coordinates
  std::vector< size_t > const & vertices = cell.entities( 0 );
  for ( size_t i = 0; i < num_vertices; ++i )
  {
    double const * coords = cell.mesh().geometry().x( vertices[i] );

    for ( size_t c = 0; c < geometric_dimension; ++c )
      coordinates[i * geometric_dimension + c] = coords[c];
  }
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_CELL_H */
