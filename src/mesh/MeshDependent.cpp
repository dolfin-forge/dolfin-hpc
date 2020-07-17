// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/MeshDependent.h>

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
MeshDependent::MeshDependent(MeshDependent const&) :
  mesh_(nullptr),
  topology_token_(0),
  geometry_token_(0)
{
}

}
