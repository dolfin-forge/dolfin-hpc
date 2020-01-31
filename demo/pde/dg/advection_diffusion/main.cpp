// Copyright (C) 2007-2008 Kristian B. Oelgaard, Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-06-29
// Last changed: 2008-07-15
//
// Steady state advection-diffusion equation, discontinuous formulation using full upwinding.

#include <dolfin.h>

#include "AdvectionDiffusion.h"
#include "OutflowFacet.h"
#include "Projection.h"

using namespace dolfin;

// Dirichlet boundary condition
struct BC : public Value<BC>
{
  void eval(real* values, const real* x) const
  {
    values[0] = sin(DOLFIN_PI*5.0*x[1]);
  }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return std::abs(x[0] - 1.0) < DOLFIN_EPS && on_boundary;
  }
};

int main(int argc, char *argv[])
{
  dolfin_init();

  // Read simple velocity field (-1.0, -0.4)
  // defined on a 64x64 unit square mesh and a quadratic vector Lagrange element
  UnitSquare mesh(64, 64);
  Function velocity(mesh);
  File("velocity.xml.gz") >> velocity;


  // Set up problem
  Matrix A;
  Vector x, b;
  Constant c(0.0); // Diffusivity constant
  Constant f(0.0); // Source term

  // Definitions for outflow facet function
  OutflowFacet::Functional M_of(mesh, velocity);

  // Penalty parameter
  Constant alpha(20.0);

  CoefficientMap bil_coef_map;
  bil_coef_map["u"] = &velocity;
  bil_coef_map["kappa"] = &c;
  bil_coef_map["alpha"] = &alpha;
  AdvectionDiffusion::BilinearForm a(mesh, bil_coef_map);
  AdvectionDiffusion::LinearForm L(mesh, f);

  // Set up boundary condition (apply strong BCs)
  Analytic<BC> g(mesh);
  DirichletBoundary boundary;
  DirichletBC bc(g, mesh, boundary, geometric);

  Assembler::assemble(A, a, true);
  Assembler::assemble(b, L, true);
  bc.apply(A, b, a);

  // Discontinuous solution
  Function uh(mesh);
  uh.init(a, 1);

  LinearSolver s;

  s.solve(A, x, b);

  // Define PDE for projection
  Projection::BilinearForm ap(mesh);
  Projection::LinearForm Lp(mesh, uh);
  // LinearPDE pde(ap, Lp, mesh);

  // Solve PDE
  // Function up(mesh);
  // pde.solve(up);

  // // Save projected solution
  // File file("temperature.pvd");
  // file << up;

  dolfin_finalize();
}
