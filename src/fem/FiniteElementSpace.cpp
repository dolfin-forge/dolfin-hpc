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
    scratch(finite_element_, dof_map_)
#if ENABLE_UFL
            ,
    ufl_(
        ufl::FiniteElementBase::create(
            ufl::Object::repr_t(finite_element_signature)))
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
    scratch(finite_element_, dof_map_)
#if ENABLE_UFL
            ,
    ufl_(ufl::FiniteElementBase::create(ufl::Object::repr_t(signature)))
#endif
{
}

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(Mesh& mesh, Form& form, uint const i) :
    mesh_(mesh),
    finite_element_(mesh.type(), form, i),
    dof_map_(DofMapCache::instance().acquire_dofmap(mesh, form, i)),
    scratch(finite_element_, dof_map_)
#if ENABLE_UFL
            ,
    ufl_(
        ufl::FiniteElementBase::create(
            ufl::Object::repr_t(element().signature())))
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
    scratch(finite_element_, dof_map_)
#if ENABLE_UFL
            ,
    ufl_(
        ufl::FiniteElementBase::create(
            ufl::Object::repr_t(element().signature())))
#endif
{
}

#if ENABLE_UFL
//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(
    Mesh& mesh, ufl::FiniteElementBase const& finite_element) :
    mesh_(mesh),
    finite_element_(finite_element),
    dof_map_(
        DofMapCache::instance().acquire_dofmap(
            mesh, DofMap::dofmap_signature(finite_element_.signature()))),
    scratch(finite_element_, dof_map_),
    ufl_(&finite_element)
{
}
#endif

//-----------------------------------------------------------------------------
FiniteElementSpace::FiniteElementSpace(FiniteElementSpace const& space,
                                       uint const& i) :
    mesh_(space.mesh()),
    finite_element_(*space.element().create_sub_element(i), true),
    dof_map_(
        DofMapCache::instance().acquire_dofmap(
            space.mesh(),
            DofMap::dofmap_signature(finite_element_.signature()))),
    scratch(finite_element_, dof_map_)
#if ENABLE_UFL
            ,
    ufl_(
        ufl::FiniteElementBase::create(
            ufl::Object::repr_t(element().signature())))
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
void FiniteElementSpace::disp() const
{
  cout << "FiniteElementSpace" << endl;
  cout << "------------------" << endl;

  // Begin indentation
  begin("");
  cout << "Finite element        : " << this->element().signature() << endl;
  cout << "Dof map               : " << this->dofmap().signature() << endl;
  // End indentation
  end();
  skip();
}

//-----------------------------------------------------------------------------
bool FiniteElementSpace::is_cellwise_defined() const
{
  return (mesh_.numCells() * dof_map_.local_dimension())
      == dof_map_.global_dimension();
}

//-----------------------------------------------------------------------------
bool FiniteElementSpace::is_cellwise_constant() const
{
  return mesh_.numCells() == dof_map_.global_dimension();
}

//-----------------------------------------------------------------------------
FiniteElementSpace::Scratch::Scratch(FiniteElement const& finite_element,
                                     DofMap const& dof_map) :
    size(0),
    space_dimension(finite_element.space_dimension()),
    local_dimension(dof_map.local_dimension()),
    dofs(NULL),
    coefficients(NULL),
    values(NULL),
    coordinates(NULL)
{
  message(1, "Creating scratch space");

  // Compute size of value (number of entries in tensor value)
  size = 1;
  for (uint i = 0; i < finite_element.value_rank(); ++i)
  {
    size *= finite_element.value_dimension(i);
  }

  // Initialize local array for mapping of dofs
  dofs = new uint[space_dimension];
  for (uint i = 0; i < space_dimension; ++i)
  {
    dofs[i] = 0;
  }

  // Initialize local array for expansion coefficients
  coefficients = new real[space_dimension];
  for (uint i = 0; i < space_dimension; ++i)
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
  coordinates = new real*[local_dimension];
  for (uint i = 0; i < local_dimension; ++i)
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
  for (uint i = 0; i < local_dimension; ++i)
  {
    delete[] coordinates[i];
  }
  delete[] coordinates;
}

//-----------------------------------------------------------------------------

}
/* namespace icorne */
