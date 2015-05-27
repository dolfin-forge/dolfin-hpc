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
    basis_(EuclideanSpace::MAX_DIMENSION, Function(mesh, 0.0)),
    node_type_(mesh)
{
}

//-----------------------------------------------------------------------------
BoundaryNormal::BoundaryNormal(BoundaryMesh& boundary) :
    mesh_(boundary.mesh()),
    boundary_(&boundary),
    local_boundary_(false),
    basis_(EuclideanSpace::MAX_DIMENSION, Function(mesh_, 0.0)),
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
void BoundaryNormal::init(FiniteElementSpace const& space)
{
  uint gdim = space.cell().dim();
  for (uint i = 0; i < EuclideanSpace::MAX_DIMENSION; ++i)
  {
    basis_[i].init(space);
  }

  if (gdim > 1)
  {
    FiniteElementSpace sub0(space, 0);
    node_type_.init(sub0);
  }
  else
  {
    node_type_.init(space);
  }
}

//-----------------------------------------------------------------------------
void BoundaryNormal::write(std::string const& filename)
{
    std::vector<std::pair<Function *, std::string> > fields;
    for(uint i = 0 ; i < EuclideanSpace::MAX_DIMENSION; ++i)
    {
      std::stringstream ss;
      ss << "E" << i;
      fields.push_back(std::pair<Function *, std::string>(&this->basis()[i],ss.str()));
    }
    fields.push_back(std::pair<Function *, std::string>(&this->node_type(),"TYPE"));
    File f(filename);
    f << fields;
}

}
