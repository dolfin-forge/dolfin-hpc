// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-02-07
// Last changed: 2007-08-20
//
// This demo program solves Poisson's equation
//
//     - div grad u(x, y) = f(x, y)
//
// on the unit square with source f given by
//
//     f(x, y) = 500*exp(-((x-0.5)^2 + (y-0.5)^2)/0.02)
//
// and boundary conditions given by
//
//     u(x, y)     = 0               for x = 0
//     du/dn(x, y) = 25 sin(5 pi y)  for x = 1
//     du/dn(x, y) = 0               otherwise

#include "Poisson.h"

#include <dolfin.h>

using namespace dolfin;

//-----------------------------------------------------------------------------

// Source term
struct Source : public Value<Source>
{
  void eval(real * value, const real* x) const
  {
    real dx = x[0] - 0.5;
    real dy = x[1] - 0.5;
    value[0] = 500.0*std::exp(-(dx*dx + dy*dy)/0.02);
  }
};

// Neumann boundary condition
struct Flux : public Value<Flux>
{
  void eval(real * value, const real* x) const
  {
    if (x[0] > DOLFIN_EPS)
      value[0] = 25.0*std::sin(5.0*DOLFIN_PI*x[1]);
    else
      value[0] = 0.0;
  }
};

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
  dolfin_init();

  // Create mesh
  Mesh mesh("UnitSquareMesh_32x32.bin");

  // Create coefficients
  Analytic<Source>  f(mesh);
  Analytic<Flux>    g(mesh);

  // Create boundary condition
  Constant u0(0.0);
  DirichletBoundary boundary;
  DirichletBC bc(u0, mesh, boundary);

  // Define PDE
  Poisson::BilinearForm a(mesh);
  Poisson::LinearForm   L(mesh, f, g);

  // Solve PDE
  Matrix A;
  Vector b;

  a.assemble(A, true);
  L.assemble(b, true);
  bc.apply(A, b, a);

  Function u(a.trial_space());

  KrylovSolver solver(bicgstab, bjacobi);
  solver.solve(A, u.vector(), b);

  message("vector l2  norm: %e", u.vector().norm());
  message("vector inf norm: %e", u.vector().max());

  // Save solution to file
  File("poisson.pvd") << u;

  return 0;
}
