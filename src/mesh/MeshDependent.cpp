// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-03-04
// Last changed: 2014-03-04

#include <dolfin/mesh/MeshDependent.h>
#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

//---------------------------------------------------------------------------
MeshDependent::MeshDependent(Mesh& mesh) :
    mesh_(mesh),
    topology_token_(0), // mesh.topology().token()
    geometry_token_(0) // mesh.geometry().token()
{
}

//---------------------------------------------------------------------------
MeshDependent::~MeshDependent()
{
}

//---------------------------------------------------------------------------
Mesh& MeshDependent::mesh() const
{
  return mesh_;
}

//---------------------------------------------------------------------------
std::string const MeshDependent::mesh_hash() const
{
  std::stringstream ss;
  ss << "Mesh@" << this << ":" << mesh_.type().description()
      << ":C" << mesh_.numCells() << ":F" << mesh_.numFacets()
      << ":V" << mesh_.numVertices() << ":T" << mesh_._timestamp;
  return ss.str();
}

//---------------------------------------------------------------------------
bool MeshDependent::invalid_mesh_topology() const
{
  return topology_token_ == mesh_.topology_token_;
}

//---------------------------------------------------------------------------
bool MeshDependent::invalid_mesh_geometry() const
{
  return geometry_token_ == mesh_.geometry_token_;
}

//---------------------------------------------------------------------------
bool MeshDependent::invalid_mesh() const
{
  return invalid_mesh_topology() || invalid_mesh_geometry();
}

//---------------------------------------------------------------------------
void MeshDependent::update_mesh_dependency()
{
  topology_token_ = mesh_.topology_token_;
  geometry_token_ = mesh_.geometry_token_;
}

}
