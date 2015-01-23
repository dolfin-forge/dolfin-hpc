// Copyright (C) 2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-11-23
// Last changed: 2008-04-28
//
// This demo program solves Poisson's equation
//
//     - div grad u(x) = f(x)
//
// on the unit interval with source f given by
//
//     f(x) = 9.0*DOLFIN_PI*DOLFIN_PI*sin(3.0*DOLFIN_PI*x[0]);
//
// and boundary conditions given by
//
//     u(x) = 0 for x = 0
//    du/dx = 0 for x = 1

#include <dolfin.h>

#include "Poisson.h"
  
using namespace dolfin;

// Boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return (std::abs(x[0]) < DOLFIN_EPS);
  }
};

// Source term
class Source : public ScalarExpression
{
public:
    
  Source() : ScalarExpression() {}

  void eval(real* values, const real* x) const
  {
    values[0] = 9.0*DOLFIN_PI*DOLFIN_PI*sin(3.0*DOLFIN_PI*x[0]);
  }

};

int main()
{
  error("Dirichlet BC gives an error in 1D");

  // Create mesh
  UnitInterval mesh(50);

  // Set up BCs
  Function zero(mesh, 0.0);
  DirichletBoundary boundary;
  DirichletBC bc(zero, mesh, boundary);

  // Create source
  Source source;
  Function f(mesh, source);

  // Define PDE
  PoissonBilinearForm a(mesh);
  PoissonLinearForm L(f);
  LinearPDE pde(a, L, mesh, bc);

  // Solve PDE
  Function u(mesh);
  pde.solve(u);

  // Save solution to file
  File file_u("poisson.pvd");
  file_u << u;

  return 0;
}
