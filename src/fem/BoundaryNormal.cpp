// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-13
// Last changed: 2014-02-13

#include <dolfin/fem/BoundaryNormal.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/mesh/BoundaryMesh.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
BoundaryNormal::BoundaryNormal(Mesh& mesh) :
    mesh_(mesh),
    boundary_(&mesh.exterior_boundary()),
    local_boundary_(false),
    basis_(3, Function(mesh)),
    node_type_(mesh)
{
}

//-----------------------------------------------------------------------------
BoundaryNormal::BoundaryNormal(BoundaryMesh& boundary) :
    mesh_(boundary.mesh()),
    boundary_(&boundary),
    local_boundary_(false),
    basis_(3, Function(mesh_)),
    node_type_(mesh_)
{
}

//-----------------------------------------------------------------------------
BoundaryNormal::~BoundaryNormal()
{
  if (local_boundary_)
  {
    delete boundary_;
  }
}

//-----------------------------------------------------------------------------
Mesh& BoundaryNormal::mesh()
{
  return mesh_;
}

//-----------------------------------------------------------------------------
BoundaryMesh& BoundaryNormal::boundary()
{
  return *boundary_;
}

//-----------------------------------------------------------------------------
Array<Function>& BoundaryNormal::basis()
{
  return basis_;
}

//-----------------------------------------------------------------------------
Function& BoundaryNormal::node_type()
{
  return node_type_;
}

//-----------------------------------------------------------------------------
void BoundaryNormal::init(Mesh& mesh, Form& form, uint i)
{
  uint gdim = mesh.geometry().dim();
  for (uint i = 0; i < 3; ++i)
  {
    basis_[i].init(mesh, form, i);
  }

  ufc::finite_element * fe = form.create_finite_element(i);
  std::string sign = fe->signature();
  if (gdim > 1)
  {
    ufc::finite_element * sub = fe->create_sub_element(0);
    sign = sub->signature();
    delete sub;
  }
  node_type_.init(mesh, sign);
  delete fe;
}

//-----------------------------------------------------------------------------
void BoundaryNormal::init(Mesh& mesh, std::string const& signature)
{
  uint gdim = mesh.geometry().dim();
  for (uint i = 0; i < 3; ++i)
  {
    basis_[i].init(mesh, signature);
  }
  if (gdim == 1)
  {
    node_type_.init(mesh, signature);
  }
  else
  {
    FiniteElement fem(signature);
    ufc::finite_element * fe = fem.create_sub_element(0);
    node_type_.init(mesh, fe->signature());
    delete fe;
  }
}

//-----------------------------------------------------------------------------
void BoundaryNormal::init(Function& other)
{
  if (other.type() != Function::discrete)
  {
    error("Initialization of BoundaryNormal from Function argument is only "
          "possible with a DiscreteFunction.");
  }
  uint gdim = other.mesh().geometry().dim();
  for (uint i = 0; i < 3; ++i)
  {
    basis_[i].init(other.mesh(), other.signature());
  }

  //TODO: Integrate this functionality to FiniteElementSpace
  if (gdim == 1)
  {
    node_type_.init(other.mesh(), other.signature());
  }
  else
  {
    ufc::finite_element * fe = other.space().element().create_sub_element(0);
    node_type_.init(other.mesh(), fe->signature());
    delete fe;
  }
}

//-----------------------------------------------------------------------------
void BoundaryNormal::write(std::string const& filename)
{
    std::vector<std::pair<Function *, std::string> > fields;
    for(uint i = 0 ; i < 3 ; ++i)
    {
      std::stringstream ss;
      ss << "E" << i;
      fields.push_back(std::pair<Function *, std::string>(this->basis()[i],ss.str()));
    }
    fields.push_back(std::pair<Function *, std::string>(this->node_type(),"TYPE"));
    File f(filename);
    f << fields;
}

}
