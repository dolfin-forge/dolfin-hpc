// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-13
// Last changed: 2014-02-13

#include <dolfin/fem/BoundaryNormal.h>
#include <dolfin/mesh/BoundaryMesh.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
BoundaryNormal::BoundaryNormal(Mesh& mesh) :
    mesh_(mesh),
    boundary_(&mesh.exterior_boundary()),
    local_boundary_(false),
    basis_(3)
{
}

//-----------------------------------------------------------------------------
BoundaryNormal::BoundaryNormal(BoundaryMesh& boundary) :
    mesh_(boundary.mesh()),
    boundary_(&boundary),
    local_boundary_(false),
    basis_(3)
{
}

//-----------------------------------------------------------------------------
BoundaryNormal::~BoundaryNormal()
{
  if(local_boundary_)
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
  uint tdim = mesh.topology().dim();
  for (uint i = 0; i < tdim; ++i)
  {
    basis_[i].init(mesh, form, i);
  }
}

//-----------------------------------------------------------------------------
void BoundaryNormal::init(Mesh& mesh, std::string const& signature)
{
  uint tdim = mesh.topology().dim();
  for (uint i = 0; i < tdim; ++i)
  {
    basis_[i].init(mesh, signature);
  }
}

//-----------------------------------------------------------------------------
void BoundaryNormal::init(Function& other)
{
  uint tdim = other.mesh().topology().dim();
  for (uint i = 0; i < tdim; ++i)
  {
    basis_[i].init(other.mesh(), other.signature());
  }
}

}
