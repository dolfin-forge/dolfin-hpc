// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson 2008.
//
// First added:  2006-06-21
// Last changed: 2008-06-26

#include <iostream>

#include <dolfin/log/log.h>
#include "BoundaryComputation.h"
#include "BoundaryMesh.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
BoundaryMesh::BoundaryMesh() : Mesh()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoundaryMesh::BoundaryMesh(Mesh& mesh) : Mesh()
{
  init(mesh);
}
//-----------------------------------------------------------------------------
BoundaryMesh::~BoundaryMesh()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void BoundaryMesh::init(Mesh& mesh)
{
  BoundaryComputation::computeBoundary(mesh, *this);
}
//-----------------------------------------------------------------------------
void BoundaryMesh::init_local(Mesh& mesh)
{
  BoundaryComputation::computeLocalBoundary(mesh, *this);
}
//-----------------------------------------------------------------------------
void BoundaryMesh::init_interior(Mesh& mesh)
{
  BoundaryComputation::computeInteriorBoundary(mesh, *this);
}
//-----------------------------------------------------------------------------

