// Copyright (C) 2008 Solveig Bruvoll and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// This demo demonstrates how to move the vertex coordinates
// of a boundary mesh and then updating the interior vertex
// coordinates of the original mesh by suitably interpolating
// the vertex coordinates (useful for implementation of ALE
// methods).

#include <dolfin.h>

using namespace dolfin;

int main()
{
  // Create mesh
  UnitSquare mesh(20, 20);

  // Create boundary mesh
  BoundaryMesh & boundary = mesh.exterior_boundary();

  // Move vertices in boundary
  for (VertexIterator v(boundary); !v.end(); ++v)
  {
    real* x = v->x();
    x[0] *= 3.0;
    x[1] += 0.1*sin(5.0*x[0]);
  }

  // Move mesh
  mesh.move(boundary);

  // Plot mesh
  File f_mesh("mesh.pvd");
  f_mesh << mesh;
  File f_mesh2("mesh.xml");
  f_mesh2 << mesh;
  File f_mesh3("boundary.xml");
  f_mesh3 << boundary;
  File f_mesh4("boundary.pvd");
  f_mesh4 << boundary;

  return 0;
}
