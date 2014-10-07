// Copyright (C) 2007-2008 Kristian B. Oelgaard, Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-06-29
// Last changed: 2008-07-15
//
// Steady state advection-diffusion equation, discontinuous formulation using full upwinding.

#include <dolfin.h>

#ifdef ENABLE_UFL
#include "ufc2/AdvectionDiffusion.h"
#include "ufc2/OutflowFacet.h"
#include "ufc2/Projection.h"
#else
#include "ufc1/AdvectionDiffusion.h"
#include "ufc1/OutflowFacet.h"
#include "ufc1/Projection.h"
#endif


using namespace dolfin;

// Dirichlet boundary condition
class BC : public ScalarExpression
{
public:

  BC() : ScalarExpression () {}

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

// Advective velocity
//class Velocity : public VectorExpression
//{
//public:
//    
//  Velocity() : VectorExpression(2) {}
//
//  void eval(real* values, const real* x) const
//  {
//    values[0] = -1.0;
//    values[1] = -0.4;
//  }
//
//  dolfin::uint rank() const
//  { return 1; }
//
//  dolfin::uint dim(dolfin::uint i) const
//  { return 2; }
//};

int main(int argc, char *argv[])
{
  // Read simple velocity field (-1.0, -0.4)
  // defined on a 64x64 unit square mesh and a quadratic vector Lagrange element
  UnitSquare mesh(64, 64);
  Function velocity(mesh);
#ifdef ENABLE_UFL 
  File vel("ufc2/velocity.xml.gz");
  vel >> velocity;
#else
  File vel("ufc1/velocity.xml.gz");
  vel >> velocity;
#endif


  // Set up problem
  Matrix A;
  Vector x, b;
  Function c(mesh, 0.0); // Diffusivity constant
  Function f(mesh, 0.0); // Source term

//  FacetNormal N(mesh);
  AvgMeshSize h(mesh);

  // Definitions for outflow facet function
  std::map<std::string const, Function *> coef_map;
  coef_map["velocity"] = &velocity;
  OutflowFacetFunctional M_of(mesh, coef_map);
  OutflowFacet of(mesh, M_of); // From SpecialFunctions.h

  // Penalty parameter
  Function alpha(mesh, 20.0);

  std::map<std::string const, Function *> bil_coef_map;
  bil_coef_map["u"] = &velocity;
  bil_coef_map["kappa"] = &c;
  bil_coef_map["alpha"] = &alpha;
  AdvectionDiffusionBilinearForm a(mesh, bil_coef_map);
  AdvectionDiffusionLinearForm L(f);

  // Set up boundary condition (apply strong BCs)
  BC g;
  Function g_func (mesh, g);
  DirichletBoundary boundary;
  DirichletBC bc(g_func, mesh, boundary, geometric);

  Assembler assembler(mesh);
  assembler.assemble(A, a, true);
  assembler.assemble(b, L, true);
  bc.apply(A, b, a);

  // Discontinuous solution
  Function uh(mesh);
  uh.init(mesh, x, a, 1);

  solve(A, x, b);

  // Define PDE for projection
  ProjectionBilinearForm ap(mesh);
  ProjectionLinearForm Lp(uh);
  LinearPDE pde(ap, Lp, mesh);

  // Solve PDE
  Function up(mesh);
  pde.solve(up);

  // Save projected solution
  File file("temperature.pvd");
  file << up;
}
