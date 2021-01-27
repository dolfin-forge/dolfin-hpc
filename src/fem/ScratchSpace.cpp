// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/ScratchSpace.h>

#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/mesh/Space.h>

#include <utility>

namespace dolfin
{

//-----------------------------------------------------------------------------

ScratchSpace::ScratchSpace( FiniteElementSpace const & space )
  : cell( space.cell() )
  , offset( 0 )
  , finite_element( &space.element() )
  , dof_map( &space.dofmap() )
  , size( finite_element->value_size() )
  , space_dimension( finite_element->space_dimension() )
  , local_dimension( dof_map->num_element_dofs() )
  , num_sub_elements( finite_element->num_sub_elements() )
  , topological_dimension( finite_element->topological_dimension() )
  , geometric_dimension( space.mesh().geometry_dimension() )
  , dofs( finite_element->space_dimension(), 0 )
  , facet_dofs( dof_map->num_facet_dofs(), 0 )
  , values( finite_element->value_size(), 0.0 )
  , coefficients( finite_element->space_dimension(), 0.0 )
  , basis_values( finite_element->space_dimension(), 0.0 )
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  , all_basis_values( finite_element->space_dimension(),
                      std::vector< real >( Space::MAX_DIMENSION, 0.0 ) )
#endif
  , coordinates( dof_map->num_element_dofs() * Space::MAX_DIMENSION, 0.0 )
  , owner_( false )
{
}

//-----------------------------------------------------------------------------

ScratchSpace::ScratchSpace( ScratchSpace const & other )
  : cell( other.cell )
  , offset( 0 )
  , finite_element( other.owner_ ? other.finite_element->create()
                                 : other.finite_element )
  , dof_map( other.owner_ ? other.dof_map->create()
                          : other.dof_map )
  , size( other.size )
  , space_dimension( other.space_dimension )
  , local_dimension( other.local_dimension )
  , num_sub_elements( other.num_sub_elements )
  , topological_dimension( other.topological_dimension )
  , geometric_dimension( other.geometric_dimension )
  , dofs( other.dofs )
  , facet_dofs( other.facet_dofs )
  , values( other.values )
  , coefficients( other.coefficients )
  , basis_values( other.basis_values )
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  , all_basis_values( other.all_basis_values )
#endif
  , coordinates( other.coordinates )
  , owner_( other.owner_ )
{
}

//-----------------------------------------------------------------------------

ScratchSpace::ScratchSpace( ScratchSpace && other )
  : cell( other.cell )
  , offset( other.offset )
  , finite_element( other.owner_ ? other.finite_element->create()
                                 : other.finite_element )
  , dof_map( other.owner_ ? other.dof_map->create()
                          : other.dof_map )
  , size( other.size )
  , space_dimension( other.space_dimension )
  , local_dimension( other.local_dimension )
  , num_sub_elements( other.num_sub_elements )
  , topological_dimension( other.topological_dimension )
  , geometric_dimension( other.geometric_dimension )
  , dofs( std::move( other.dofs ) )
  , facet_dofs( std::move( other.facet_dofs ) )
  , values( std::move( other.values ) )
  , coefficients( std::move( other.coefficients ) )
  , basis_values( std::move( other.basis_values ) )
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  , all_basis_values( std::move( other.all_basis_values ) )
#endif
  , coordinates( std::move( other.coordinates ) )
  , owner_( other.owner_ )
{
}

//-----------------------------------------------------------------------------

ScratchSpace::ScratchSpace( FiniteElementSpace const & space,
                            SubSystem const &          sub_system )
  : cell( space.cell() )
  , offset( 0 )
  , finite_element( space.element().create_sub_element( sub_system.array() ) )
  , dof_map( space.dofmap().create_sub_dofmap( sub_system.array(), offset ) )
  , size( finite_element->value_size() )
  , space_dimension( finite_element->space_dimension() )
  , local_dimension( dof_map->num_element_dofs() )
  , num_sub_elements( finite_element->num_sub_elements() )
  , topological_dimension( finite_element->topological_dimension() )
  , geometric_dimension( space.mesh().geometry_dimension() )
  , dofs( finite_element->space_dimension(), 0 )
  , facet_dofs( dof_map->num_facet_dofs(), 0 )
  , values( finite_element->value_size(), 0.0 )
  , coefficients( finite_element->space_dimension(), 0.0 )
  , basis_values( finite_element->space_dimension(), 0.0 )
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  , all_basis_values( finite_element->space_dimension(),
                      std::vector< real >( Space::MAX_DIMENSION, 0.0 ) )
#endif
  , coordinates( dof_map->num_element_dofs() * Space::MAX_DIMENSION, 0.0 )
  , owner_( true )
{
}

//-----------------------------------------------------------------------------

ScratchSpace::~ScratchSpace()
{
  if ( owner_ )
  {
    delete dof_map;
    delete finite_element;
  }
}

//-----------------------------------------------------------------------------

ScratchSpace & ScratchSpace::operator=( ScratchSpace & other )
{
  dolfin_assert( other.finite_element->signature() == finite_element->signature() );
  dolfin_assert( other.dof_map->signature()        == dof_map->signature() );

  cell         = other.cell;
  offset       = other.offset;
  dofs         = other.dofs;
  facet_dofs   = other.facet_dofs;
  values       = other.values;
  coefficients = other.coefficients;
  basis_values = other.basis_values;
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  all_basis_values = other.all_basis_values;
#endif
  coordinates  = other.coordinates;

  return *this;
}

//-----------------------------------------------------------------------------

ScratchSpace & ScratchSpace::operator=( ScratchSpace && other )
{
  dolfin_assert( other.finite_element->signature() == finite_element->signature() );
  dolfin_assert( other.dof_map->signature()        == dof_map->signature() );

  cell         = std::move( other.cell );
  offset       = std::move( other.offset );
  dofs         = std::move( other.dofs );
  facet_dofs   = std::move( other.facet_dofs );
  values       = std::move( other.values );
  coefficients = std::move( other.coefficients );
  basis_values = std::move( other.basis_values );
#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  all_basis_values = std::move( other.all_basis_values =;
#endif
  coordinates  = std::move( other.coordinates );

  return *this;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
