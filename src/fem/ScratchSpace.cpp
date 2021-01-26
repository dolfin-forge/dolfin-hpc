// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/ScratchSpace.h>

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/mesh/Space.h>
#include <dolfin/mesh/entities/Cell.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

ScratchSpace::ScratchSpace( FiniteElementSpace const & space )
  : cell( space.cell() )
  , offset( 0 )
  , finite_element( &space.element() )
  , dof_map( &space.dofmap() )
  , size( value_size( *finite_element ) )
  , space_dimension( finite_element->space_dimension() )
  , local_dimension( dof_map->num_element_dofs() )
  , num_sub_elements( finite_element->num_sub_elements() )
  , topological_dimension( finite_element->topological_dimension() )
  , geometric_dimension( space.mesh().geometry_dimension() )
  , dofs( new size_t[space_dimension] )
  , facet_dofs( new size_t[dof_map->num_facet_dofs()] )
  , values( new real[size] )
  , coefficients( new real[space_dimension] )
  , basis_values( new real[space_dimension] )
  ,
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  all_basis_values( new real *[space_dimension] )
  ,
#endif
  coordinates( dof_map->num_element_dofs() * Space::MAX_DIMENSION, 0.0 )
  , owner_( false )
{
  init();
}

//-----------------------------------------------------------------------------

ScratchSpace::ScratchSpace( FiniteElementSpace const & space,
                            SubSystem const &          sub_system )
  : cell( space.cell() )
  , offset( 0 )
  , finite_element( space.element().create_sub_element( sub_system.array() ) )
  , dof_map( space.dofmap().create_sub_dofmap( sub_system.array(), offset ) )
  , size( value_size( *finite_element ) )
  , space_dimension( finite_element->space_dimension() )
  , local_dimension( dof_map->num_element_dofs() )
  , num_sub_elements( finite_element->num_sub_elements() )
  , topological_dimension( finite_element->topological_dimension() )
  , geometric_dimension( space.mesh().geometry_dimension() )
  , dofs( new size_t[space_dimension] )
  , facet_dofs( new size_t[dof_map->num_facet_dofs()] )
  , values( new real[size] )
  , coefficients( new real[space_dimension] )
  , basis_values( new real[space_dimension] )
  ,
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  all_basis_values( new real *[space_dimension] )
  ,
#endif
  coordinates( dof_map->num_element_dofs() * Space::MAX_DIMENSION, 0.0 )
  , owner_( true )
{
  init();
}

//-----------------------------------------------------------------------------

ScratchSpace::ScratchSpace( ScratchSpace const & other )
  : cell( other.cell )
  , offset( 0 )
  , finite_element( nullptr )
  , dof_map( nullptr )
  , size( 0 )
  , space_dimension( 0 )
  , local_dimension( 0 )
  , num_sub_elements( 0 )
  , topological_dimension( 0 )
  , geometric_dimension( 0 )
  , dofs( nullptr )
  , facet_dofs( nullptr )
  , values( nullptr )
  , coefficients( nullptr )
  , basis_values( nullptr )
  ,
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  all_basis_values( nullptr )
  ,
#endif
  coordinates( 0 )
  , owner_( false )
{
  error( "ScratchSpace::ScratchSpace(ScratchSpace const& other)" );
}

//-----------------------------------------------------------------------------

ScratchSpace::~ScratchSpace()
{
  coordinates.clear();
  delete[] basis_values;
  delete[] coefficients;
  delete[] values;
  delete[] facet_dofs;
  delete[] dofs;
  if ( owner_ )
  {
    delete dof_map;
    delete finite_element;
  }
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  for ( size_t i = 0; i < space_dimension; ++i )
  {
    if ( all_basis_values[i] != nullptr )
      delete[] all_basis_values[i];
  }
  delete[] all_basis_values;
#endif
}

//-----------------------------------------------------------------------------

auto ScratchSpace::value_size( ufc::finite_element const & finite_element )
  -> size_t
{
  // Compute size of value (number of entries in tensor value)
  size_t size = 1;
  for ( size_t i = 0; i < finite_element.value_rank(); ++i )
  {
    size *= finite_element.value_dimension( i );
  }
  return size;
}

//-----------------------------------------------------------------------------

void ScratchSpace::init()
{
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  // Initialize local array for dof coordinates
  for ( size_t i = 0; i < space_dimension; ++i )
  {
    // Using same storage size as a Point
    all_basis_values[i] = new real[Space::MAX_DIMENSION];
    std::fill_n( all_basis_values[i], Space::MAX_DIMENSION, 0.0 );
  }
#endif
}

} /* namespace dolfin */
