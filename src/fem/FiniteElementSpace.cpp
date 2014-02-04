// Copyright (C) 2013 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-03
// Last changed: 2014-02-03

#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/DofMapCache.h>

#include <dolfin/fem/DofMap.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(
    Mesh& mesh, std::string const& finite_element_signature,
    std::string const& dof_map_signature) :
    mesh_(mesh),
    finite_element_(finite_element_signature),
    dof_map_(DofMapCache::instance().acquire_dofmap(mesh, dof_map_signature)),
    scratch(finite_element_)
#if ENABLE_UFL
            ,
    ufl_class_(ufl::Object::repr_t(element().signature()))
#endif
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Mesh& mesh, std::string const& signature) :
    mesh_(mesh),
    finite_element_(signature),
    dof_map_(
        DofMapCache::instance().acquire_dofmap(
            mesh, DofMap::dofmap_signature(signature))),
    scratch(finite_element_)
#if ENABLE_UFL
            ,
    ufl_class_(ufl::Object::repr_t(element().signature()))
#endif
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Mesh& mesh, Form& form, uint const i) :
    mesh_(mesh),
    finite_element_(mesh.type(), form, i),
    dof_map_(DofMapCache::instance().acquire_dofmap(mesh, form, i)),
    scratch(finite_element_)
#if ENABLE_UFL
            ,
    ufl_class_(ufl::Object::repr_t(element().signature()))
#endif
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Mesh& mesh,
                                       ufc::finite_element& finite_element,
                                       bool const finite_element_local) :
    mesh_(mesh),
    finite_element_(finite_element, finite_element_local),
    dof_map_(
        DofMapCache::instance().acquire_dofmap(
            mesh, DofMap::dofmap_signature(finite_element_.signature()))),
    scratch(finite_element_)
#if ENABLE_UFL
            ,
    ufl_class_(ufl::Object::repr_t(element().signature()))
#endif
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(FiniteElementSpace const& space,
                                       uint const& i) :
    mesh_(space.mesh()),
    finite_element_(*space.element().create_sub_element(i), true),
    dof_map_(
        DofMapCache::instance().acquire_dofmap(
            space.mesh(),
            DofMap::dofmap_signature(finite_element_.signature()))),
    scratch(finite_element_)
#if ENABLE_UFL
            ,
    ufl_class_(ufl::Object::repr_t(element().signature()))
#endif
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::~FiniteElementSpace()
{
  DofMapCache::instance().release_dofmap(dof_map_);
}

//-----------------------------------------------------------------------------
FiniteElement const& FiniteElementSpace::element() const
{
  return finite_element_;
}

//-----------------------------------------------------------------------------
DofMap const& FiniteElementSpace::dofmap() const
{
  return dof_map_;
}

//-----------------------------------------------------------------------------
FiniteElementSpace::Scratch::Scratch(FiniteElement const& finite_element) :
    size(0),
    dimension(finite_element.space_dimension()),
    dofs(NULL),
    coefficients(NULL),
    values(NULL),
    coordinates(NULL)
{
  // Compute size of value (number of entries in tensor value)
  size = 1;
  for (uint i = 0; i < finite_element.value_rank(); ++i)
  {
    size *= finite_element.value_dimension(i);
  }

  // Initialize local array for mapping of dofs
  dofs = new uint[dimension];
  for (uint i = 0; i < dimension; ++i)
  {
    dofs[i] = 0;
  }

  // Initialize local array for expansion coefficients
  coefficients = new real[dimension];
  for (uint i = 0; i < dimension; ++i)
  {
    coefficients[i] = 0.0;
  }

  // Initialize local array for values
  values = new real[size];
  for (uint i = 0; i < size; ++i)
  {
    values[i] = 0.0;
  }

  // Initialize local array for dof coordinates
  for (uint i = 0; i < dimension; ++i)
  {
    coordinates[i] = new real[3]; // Internally Point is implemented for d = 3
  }

}

//-----------------------------------------------------------------------------
FiniteElementSpace::Scratch::~Scratch()
{
  delete[] dofs;
  delete[] coefficients;
  delete[] values;
  for (uint i = 0; i < dimension; ++i)
  {
    delete[] coordinates[i];
  }
  delete[] coordinates;
}

//-----------------------------------------------------------------------------

}
/* namespace icorne */
