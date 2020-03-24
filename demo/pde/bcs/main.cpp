// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-05-23
// Last changed: 2008-05-23
//
// This demo illustrates how to set boundary conditions for meshes
// that include boundary indicators. The mesh used in this demo was
// generated with VMTK (http://villacamozzi.marionegri.it/~luca/vmtk/).

#include <dolfin.h>
#include "ufc1/Poisson.h"

using namespace dolfin;

//-----------------------------------------------------------------------------

// Sub domain for Dirichlet boundary condition
struct DirichletBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS && on_boundary;
  }
};

//-----------------------------------------------------------------------------


int main()
{
  // Create mesh and finite element
  Mesh mesh("../../../data/meshes/aneurysm.bin");

  // Define variational problem
  Constant f(0.0);
  Poisson::BilinearForm a(mesh);
  Poisson::LinearForm L(mesh,f);

  // Define boundary condition values
  Constant u0(0.0);
  Constant u1(1.0);
  Constant u2(2.0);
  Constant u3(3.0);

  // Define boundary conditions
  DirichletBoundary boundary;
  DirichletBC bc0(u0, mesh, boundary);
  DirichletBC bc1(u1, mesh, boundary);
  DirichletBC bc2(u2, mesh, boundary);
  DirichletBC bc3(u3, mesh, boundary);

  // Solve PDE
  Matrix A;
  Vector b;

  a.assemble( A, true );
  L.assemble( b, true );
  bc0.apply( A, b, a );
  bc1.apply( A, b, a );
  bc2.apply( A, b, a );
  bc3.apply( A, b, a );

  Function     u( a.trial_space() );
  LinearSolver solver;

  solver.solve( A, u.vector(), b );
  u.sync();

  // Save solution to VTK format
  File vtk_file( "bcs.pvd" );
  vtk_file << u;

  return 0;
}
