// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007, 2008.
// Modified by Aurélien Larcher, 2013.
//
// First added:  2008-06-18
// Last changed: 2013-09-13

#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/fem/BoundaryCondition.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
BoundaryCondition::BoundaryCondition(std::string const& type,
                                     Mesh& mesh,
                                     SubDomain const& sub_domain) :
    type_(type),
    mesh_(mesh),
    data_(NULL),
    sub_domain_index_(0),
    has_geometrical_sub_domain_(true),
    geometrical_sub_domain_(&sub_domain),
    sub_domain_markers_(),
    local_sub_domain_markers_(true),
    sub_system_()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoundaryCondition::BoundaryCondition(std::string const& type,
                                     MeshFunction<uint>& sub_domains,
                                     uint sub_domain) :
    type_(type),
    mesh_(sub_domains.mesh()),
    data_(NULL),
    sub_domain_index_(sub_domain),
    has_geometrical_sub_domain_(false),
    geometrical_sub_domain_(NULL),
    sub_domain_markers_(&sub_domains),
    local_sub_domain_markers_(false),
    sub_system_()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoundaryCondition::BoundaryCondition(std::string const& type,
                                     Mesh& mesh,
                                     SubDomain const& sub_domain,
                                     SubSystem const sub_system) :
    type_(type),
    mesh_(mesh),
    data_(NULL),
    sub_domain_index_(0),
    has_geometrical_sub_domain_(true),
    geometrical_sub_domain_(&sub_domain),
    sub_domain_markers_(),
    local_sub_domain_markers_(true),
    sub_system_(sub_system)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoundaryCondition::BoundaryCondition(std::string const& type,
                                     MeshFunction<uint>& sub_domains,
                                     uint sub_domain,
                                     SubSystem const sub_system) :
    type_(type),
    mesh_(sub_domains.mesh()),
    data_(NULL),
    sub_domain_index_(sub_domain),
    has_geometrical_sub_domain_(false),
    geometrical_sub_domain_(NULL),
    sub_domain_markers_(&sub_domains),
    local_sub_domain_markers_(false),
    sub_system_(sub_system)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoundaryCondition::~BoundaryCondition()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoundaryCondition::LocalData::LocalData(Mesh& mesh, SubSystem const& sub_system,
                                        Form const& form) :
    ufc_mesh(mesh),
    finite_element(NULL),
    dof_map(NULL),
    offset(0),
    w(NULL),
    cell_dofs(NULL),
    facet_dofs(NULL),
    is_subspace_(sub_system.depth() > 0)
{
  // Check arity of form
  if (form.rank() != 2)
  {
    error("Form must be bilinear for application of boundary conditions.");
  }

  // Create finite element (second argument of form)
  finite_element = form.create_finite_element(1);

  // Extract sub element and sub dof map if we have a sub system
  if (is_subspace_)
  {
    ufc::finite_element * sub_finite_element =
      FiniteElement::create_sub_element(*finite_element, sub_system.array());
    delete finite_element;
    finite_element = sub_finite_element;

    dof_map = new DofMap(form.dofmaps()[1], sub_system.array(), offset);
  }
  else
  {
    dof_map = &form.dofmaps()[1];
  }

  // Create local data used to set boundary conditions
  w = new real[finite_element->space_dimension()];
  cell_dofs = new uint[finite_element->space_dimension()];
  for (uint i = 0; i < finite_element->space_dimension(); i++)
  {
    w[i] = 0.0;
    cell_dofs[i] = 0;
  }
  facet_dofs = new uint[dof_map->num_facet_dofs()];
  for (uint i = 0; i < dof_map->num_facet_dofs(); ++i)
  {
    facet_dofs[i] = 0;
  }

  // Create local coordinate data
  coordinates = new real*[dof_map->local_dimension()];
  for (uint i = 0; i < dof_map->local_dimension(); ++i)
  {
    coordinates[i] = new real[mesh.geometry().dim()];
    for (uint j = 0; j < mesh.geometry().dim(); ++j)
    {
      coordinates[i][j] = 0.0;
    }
  }
}
//-----------------------------------------------------------------------------
BoundaryCondition::LocalData::~LocalData()
{
  if (coordinates)
  {
    for (uint i = 0; i < dof_map->local_dimension(); ++i)
    {
      delete[] coordinates[i];
    }
    delete[] coordinates;
  }

  //TODO: Always delete it until we fix the problem of the unkown space at the
  //      Form level
  delete finite_element;

  if(is_subspace_)
  {
    delete dof_map;
  }

  delete[] w;
  delete[] cell_dofs;
  delete[] facet_dofs;
}
//-----------------------------------------------------------------------------

}
