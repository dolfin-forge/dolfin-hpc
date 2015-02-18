// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-07-11
// Last changed: 2007-08-20
//
// This demo program solves Poisson's equation
//
//     - div grad u(x, y) = f(x, y)
//
// on the unit square with homogeneous Dirichlet boundary conditions
// at y = 0, 1 and periodic boundary conditions at x = 0, 1.

#include <dolfin.h>
#include "Poisson.h"
  
using namespace dolfin;

int main()
{
  // Source term
  class Source : public Expression
  {
  public:
    
    Source() : Expression() {}

    void eval(real* values, const real* x) const
    {
      real dx = x[0] - 0.5;
      real dy = x[1] - 0.5;
      values[0] = x[0]*sin(5.0*DOLFIN_PI*x[1]) + 1.0*exp(-(dx*dx + dy*dy)/0.02);
    }

    uint dim(uint i) const
    {
      return 1;
    }

    uint rank() const
    {
      return 0;
    }
  };

  // Sub domain for Dirichlet boundary condition
  class DirichletBoundary : public SubDomain
  {
    bool inside(const real* x, bool on_boundary) const
    {
      return (x[1] < DOLFIN_EPS || x[1] > (1.0 - DOLFIN_EPS)) && on_boundary;
    }
  };

  // Sub domain for Periodic boundary condition
  class PeriodicBoundary : public SubDomain
  {
    bool inside(const real* x, bool on_boundary) const
    {
      return x[0] < DOLFIN_EPS && x[0] > -DOLFIN_EPS && on_boundary;
    }

    void map(const real* x, real* y) const
    {
      y[0] = x[0] - 1.0;
      y[1] = x[1];
    }
  };

  // Create mesh
  UnitSquare mesh(32, 32);

  // Create functions
  Source source;
  Function f(mesh, source);

  // Create Dirichlet boundary condition
  Function u0(mesh, 0.0);
  DirichletBoundary dirichlet_boundary;
  DirichletBC bc0(u0, mesh, dirichlet_boundary);
  
  // Create periodic boundary condition
  PeriodicBoundary periodic_boundary;
  PeriodicBC bc1(mesh, periodic_boundary);

  // Collect boundary conditions
  Array<BoundaryCondition*> bcs;
  bcs.push_back(&bc0);
  bcs.push_back(&bc1);

  // Define PDE
  PoissonBilinearForm a(mesh);
  PoissonLinearForm L(f);

  LinearPDE pde(a, L, mesh, bcs);

  // Solve PDE
  Function u;
  pde.solve(u);

  // Save solution to file
  File file("poisson.pvd");
  file << u;

  return 0;
}
