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
#include <dolfin/mesh/PeriodicSubDomain.h>

using namespace dolfin;

// Source term
struct Source : public Value< Source >
{
  void eval(real* values, const real* x) const
  {
    real dx = x[0] - 0.5;
    real dy = x[1] - 0.5;
    values[0] = x[0]*sin(5.0*DOLFIN_PI*x[1]) + 1.0*exp(-(dx*dx + dy*dy)/0.02);
  }
};

// Sub domain for Dirichlet boundary condition
struct DirichletBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return (x[1] < DOLFIN_EPS || x[1] > (1.0 - DOLFIN_EPS)) && on_boundary;
  }
};

// Sub domain for Periodic boundary condition
struct PeriodicBoundary : public PeriodicSubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS && x[0] > -DOLFIN_EPS && on_boundary;
  }

  void map(const real* xH, real* xG) const
  {
    xG[0] = xH[0] - 1.0;
    xG[1] = xH[1];
  }
};

class PeriodicBC1 : public PeriodicBC
{
public:
  PeriodicBC1(Mesh& mesh, PeriodicSubDomain const& sub_domain)
    : PeriodicBC(mesh, sub_domain)
    {}


  inline void sync(Time const& t) {}
};

int main()
{
  // Create mesh
  UnitSquare mesh(32, 32);

  // Create functions
  Analytic<Source>  f(mesh);

  // Create Dirichlet boundary condition
  Constant u0(0.0);
  DirichletBoundary dirichlet_boundary;
  DirichletBC bc0(u0, mesh, dirichlet_boundary);

  // Create periodic boundary condition
  PeriodicBoundary periodic_boundary;
  PeriodicBC1 bc1(mesh, periodic_boundary);

  // Define PDE
  Poisson::BilinearForm a(mesh);
  Poisson::LinearForm L(mesh,f);

  // Solve PDE
  Matrix A;
  Vector b;

  a.assemble(A, true);
  L.assemble(b, true);
  bc0.apply(A, b, a);
  bc1.apply(A, b, a);


  Function u(a.trial_space());
  KrylovSolver solver(bicgstab, bjacobi);

  solver.solve(A, u.vector(), b);
  u.sync();

  // Save solution to file
  File file("poisson.pvd");
  file << u;

  return 0;
}
