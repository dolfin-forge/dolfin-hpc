// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-13
// Last changed: 2014-02-13

#include <dolfin/fem/ScratchSpace.h>

#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Space.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
ScratchSpace::ScratchSpace(FiniteElementSpace const& space) :
    mesh(space.mesh()),
    cell(space.cell()),
    offset(0),
    finite_element(&space.element()),
    dof_map(&space.dofmap()),
    size(value_size(*finite_element)),
    space_dimension(finite_element->space_dimension()),
    local_dimension(dof_map->local_dimension()),
    num_sub_elements(finite_element->num_sub_elements()),
    topological_dimension(finite_element->topological_dimension()),
    geometric_dimension(space.mesh().geometry().dim()),
    dofs(new uint[space_dimension]),
    facet_dofs(new uint[dof_map->num_facet_dofs()]),
    values(new real[size]),
    coefficients(new real[space_dimension]),
    basis_values(new real[space_dimension]),
    coordinates(new real*[local_dimension]),
    owner_(false)
{
  init();
}

//-----------------------------------------------------------------------------
ScratchSpace::ScratchSpace(FiniteElementSpace const& space,
                           SubSystem const& sub_system) :
    mesh(space.mesh()),
    cell(space.cell()),
    offset(0),
    finite_element(space.element().create_sub_element(sub_system.array())),
    dof_map(space.dofmap().create_sub_dofmap(sub_system.array(), offset)),
    size(value_size(*finite_element)),
    space_dimension(finite_element->space_dimension()),
    local_dimension(dof_map->local_dimension()),
    num_sub_elements(finite_element->num_sub_elements()),
    topological_dimension(finite_element->topological_dimension()),
    geometric_dimension(space.mesh().geometry().dim()),
    dofs(new uint[space_dimension]),
    facet_dofs(new uint[dof_map->num_facet_dofs()]),
    values(new real[size]),
    coefficients(new real[space_dimension]),
    basis_values(new real[space_dimension]),
    coordinates(new real*[local_dimension]),
    owner_(true)
{
  init();
}

//-----------------------------------------------------------------------------
ScratchSpace::ScratchSpace(ScratchSpace const& other) :
    mesh(other.mesh),
    cell(other.cell),
    offset(0),
    finite_element(NULL),
    dof_map(NULL),
    size(0),
    space_dimension(0),
    local_dimension(0),
    num_sub_elements(0),
    topological_dimension(0),
    geometric_dimension(0),
    dofs(NULL),
    facet_dofs(NULL),
    values(NULL),
    coefficients(NULL),
    basis_values(NULL),
    coordinates(NULL),
    owner_(false)
{
  error("ScratchSpace::ScratchSpace(ScratchSpace const& other)");
}

//-----------------------------------------------------------------------------
ScratchSpace::~ScratchSpace()
{
  for (uint i = 0; i < local_dimension; ++i)
  {
    delete[] coordinates[i];
  }
  delete[] coordinates;
  delete[] basis_values;
  delete[] coefficients;
  delete[] values;
  delete[] facet_dofs;
  delete[] dofs;
  if (owner_)
  {
    delete dof_map;
    delete finite_element;
  }
}

//-----------------------------------------------------------------------------
uint ScratchSpace::value_size(ufc::finite_element const& finite_element)
{
  // Compute size of value (number of entries in tensor value)
  uint size = 1;
  for (uint i = 0; i < finite_element.value_rank(); ++i)
  {
    size *= finite_element.value_dimension(i);
  }
  return size;
}

//-----------------------------------------------------------------------------
void ScratchSpace::init()
{
  // Initialize local array for dof coordinates
  for (uint i = 0; i < local_dimension; ++i)
  {
    // Using same storage size as a Point
    coordinates[i] = new real[Space::MAX_DIMENSION];
    std::fill_n(coordinates[i],  Space::MAX_DIMENSION, 0.0);
  }
}

} /* namespace dolfin */
