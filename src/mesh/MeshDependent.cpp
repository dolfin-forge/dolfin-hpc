// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/MeshDependent.h>
#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

//---------------------------------------------------------------------------
MeshDependent::MeshDependent(Mesh& mesh) :
  mesh_(&mesh),
  topology_token_(mesh.topology().token()),
  geometry_token_(mesh.geometry().token())
{
}

//---------------------------------------------------------------------------
MeshDependent::MeshDependent(MeshDependent const& other) :
  mesh_(NULL),
  topology_token_(0),
  geometry_token_(0)
{
}

//---------------------------------------------------------------------------
MeshDependent::~MeshDependent()
{
}

//---------------------------------------------------------------------------
Mesh& MeshDependent::mesh() const
{
  return *mesh_;
}

//---------------------------------------------------------------------------
bool MeshDependent::invalid_mesh_topology() const
{
  return topology_token_ != mesh_->topology().token();
}

//---------------------------------------------------------------------------
bool MeshDependent::invalid_mesh_geometry() const
{
  return geometry_token_ != mesh_->geometry().token();
}

//---------------------------------------------------------------------------
bool MeshDependent::invalid_mesh() const
{
  return invalid_mesh_topology() || invalid_mesh_geometry();
}

//---------------------------------------------------------------------------
void MeshDependent::update_mesh_dependency()
{
  topology_token_ = mesh_->topology().token();
  geometry_token_ = mesh_->geometry().token();
}

}
