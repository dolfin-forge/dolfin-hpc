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
#include <dolfin/mesh/EuclideanSpace.h>

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
    dofs(new uint[space_dimension]),
    facet_dofs(new uint[dof_map->num_facet_dofs()]),
    coefficients(new real[space_dimension]),
    values(new real[size]),
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
    dofs(new uint[space_dimension]),
    facet_dofs(new uint[dof_map->num_facet_dofs()]),
    coefficients(new real[space_dimension]),
    values(new real[size]),
    coordinates(new real*[local_dimension]),
    owner_(true)
{
  init();
}

//-----------------------------------------------------------------------------
ScratchSpace::~ScratchSpace()
{
  for (uint i = 0; i < local_dimension; ++i)
  {
    delete[] coordinates[i];
  }
  delete[] coordinates;
  delete[] values;
  delete[] coefficients;
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
  message(1, "Creating scratch space");

  // Initialize local array for dof coordinates
  for (uint i = 0; i < local_dimension; ++i)
  {
    // Using same storage size as a Point
    coordinates[i] = new real[EuclideanSpace::MAX_DIMENSION];
    std::memset(&coordinates[i][0], 0.0,
                EuclideanSpace::MAX_DIMENSION * sizeof(real));
  }
}

} /* namespace dolfin */
