// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-13
// Last changed: 2014-02-13

#include <dolfin/fem/ScratchSpace.h>

#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/mesh/Cell.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
ScratchSpace::ScratchSpace(Cell const& cell,
                           FiniteElement const& finite_element,
                           DofMap const& dof_map) :
    size(value_size(finite_element)),
    space_dimension(finite_element.space_dimension()),
    local_dimension(dof_map.local_dimension()),
    num_sub_elements(finite_element.num_sub_elements()),
    topological_dimension(finite_element.topological_dimension()),
    dofs(new uint[space_dimension]),
    coefficients(new real[space_dimension]),
    values(new real[size]),
    coordinates(new real*[local_dimension]),
    tabulation_on_cell_(new uint[space_dimension]),
    tabulation_per_sub_element_(new uint*[num_sub_elements]),
    sub_element_space_dimensions_(new uint[num_sub_elements]),
    tabulation_per_entity_(new uint*[topological_dimension + 1]),
    num_entity_dofs_(new uint[topological_dimension + 1])
{
  Initialize(cell, finite_element, dof_map);
}

//-----------------------------------------------------------------------------
ScratchSpace::ScratchSpace(FiniteElementSpace const& space) :
    size(value_size(space.element())),
    space_dimension(space.element().space_dimension()),
    local_dimension(space.dofmap().local_dimension()),
    num_sub_elements(space.element().num_sub_elements()),
    topological_dimension(space.element().topological_dimension()),
    dofs(new uint[space_dimension]),
    coefficients(new real[space_dimension]),
    values(new real[size]),
    coordinates(new real*[local_dimension]),
    tabulation_on_cell_(new uint[space_dimension]),
    tabulation_per_sub_element_(new uint*[num_sub_elements]),
    sub_element_space_dimensions_(new uint[num_sub_elements]),
    tabulation_per_entity_(new uint*[topological_dimension + 1]),
    num_entity_dofs_(new uint[topological_dimension + 1])
{
  Initialize(space.cell(), space.element(), space.dofmap());
}

//-----------------------------------------------------------------------------
ScratchSpace::~ScratchSpace()
{
  delete[] num_entity_dofs_;
  for (uint i = 0; i < topological_dimension + 1; ++i)
  {
    delete[] tabulation_per_entity_[i];
  }
  delete[] tabulation_per_entity_;
  delete[] sub_element_space_dimensions_;
  for (uint i = 0; i < num_sub_elements; ++i)
  {
    delete[] tabulation_per_sub_element_[i];
  }
  delete[] tabulation_per_sub_element_;
  delete[] tabulation_on_cell_;
  for (uint i = 0; i < local_dimension; ++i)
  {
    delete[] coordinates[i];
  }
  delete[] coordinates;
  delete[] values;
  delete[] coefficients;
  delete[] dofs;
}

//-----------------------------------------------------------------------------
void ScratchSpace::Initialize(Cell const& cell,
                              FiniteElement const& finite_element,
                              DofMap const& dof_map)
{
  message(1, "Creating scratch space");

  // Initialize local array for dof coordinates
  for (uint i = 0; i < local_dimension; ++i)
  {
    coordinates[i] = new real[3]; // Internally Point is implemented for d = 3
  }

  // Initialize tabulation per subspace
  for (uint s = 0; s < num_sub_elements; ++s)
  {
    ufc::finite_element * sub = finite_element.create_sub_element(s);
    sub_element_space_dimensions_[s] = sub->space_dimension();
    tabulation_per_sub_element_[s] = new uint[sub_element_space_dimensions_[s]];
    delete sub;
  }

  // Initialize tabulation per entity
  for (uint d = 0; d < topological_dimension + 1; ++d)
  {
    num_entity_dofs_[d] = dof_map.num_entity_dofs(d);
    tabulation_per_entity_[d] = new uint[num_entity_dofs_[d]];
  }
}

//-----------------------------------------------------------------------------
uint ScratchSpace::value_size(FiniteElement const& finite_element)
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
void ScratchSpace::set_cell_tabulation(Cell const& cell,
                                       ufc::dof_map const& dof_map,
                                       uint **& dofs)
{
  for (uint e = 0; e < topological_dimension; ++e)
  {
    uint *& entity_dofs = dofs[e];
    for()
    {
      tabulate_entity_dofs()
    }
  }
}

} /* namespace dolfin */
