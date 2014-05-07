// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson 2008.
//
// First added:  2006-06-21
// Last changed: 2008-06-26

#include <iostream>

#include <dolfin/log/log.h>
#include <dolfin/mesh/BoundaryComputation.h>
#include <dolfin/mesh/BoundaryMesh.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
BoundaryMesh::BoundaryMesh(Mesh& mesh, BoundaryMesh::Type type) :
    Mesh(),
    mesh_(mesh),
    global_mesh_hash_(mesh.hash())
{
  switch (type)
  {
  case BoundaryMesh::exterior:
    // Exterior boundary i.e facets at the domain boundary
    BoundaryComputation::computeBoundary(mesh, *this);
    break;
  case BoundaryMesh::interior:
    // Interior boundary i.e facets between processors
    BoundaryComputation::computeInteriorBoundary(mesh, *this);
    break;
  case BoundaryMesh::full:
    // Full boundary incl. facets between processors
    BoundaryComputation::computeLocalBoundary(mesh, *this);
    break;
  default:
    error("Unknown boundary mesh type.");
    break;
  }
}

//-----------------------------------------------------------------------------
BoundaryMesh::~BoundaryMesh()
{
  // Do nothing
}

//-----------------------------------------------------------------------------
Mesh& BoundaryMesh::mesh()
{
  return mesh_;
}

//-----------------------------------------------------------------------------
std::string const BoundaryMesh::mesh_hash()
{
  return global_mesh_hash_;
}

}
