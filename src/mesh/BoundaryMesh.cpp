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
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
BoundaryMesh::BoundaryMesh(Mesh& mesh, BoundaryMesh::Type type) :
    Mesh(),
    MeshDependent(mesh),
    type_(type)
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
uint BoundaryMesh::facet_index(Cell const& boundary_cell)
{
  return this->data().meshFunction("cell map")->get(boundary_cell);
}

//-----------------------------------------------------------------------------
uint BoundaryMesh::vertex_index(Vertex const& boundary_vertex)
{
  return this->data().meshFunction("vertex map")->get(boundary_vertex);
}

//-----------------------------------------------------------------------------
BoundaryMesh::Type BoundaryMesh::boundary_type() const
{
  return type_;
}

//-----------------------------------------------------------------------------

}
