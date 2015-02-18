// Copyright (C) 2006-2007 Anders Logg and Kristian Oelgaard.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-12-05
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
//     u(x, y)     = 0
//     du/dn(x, y) = 0
//
// using a discontinuous Galerkin formulation (interior penalty method).

#include <dolfin.h>

#include "Poisson.h"
#include "P1Projection.h"

using namespace dolfin;

int main()
{
  // Source term
  class Source : public ScalarExpression
  {
  public:
    
    Source() : ScalarExpression() {}

    void eval(real* values, const real* x) const
    {
      real dx = x[0] - 0.5;
      real dy = x[1] - 0.5;
      values[0] = 500.0*exp(-(dx*dx + dy*dy)/0.02);
    }
  };
 
  // Create mesh
  UnitSquare mesh(64, 64);

  // Create functions
  Source f;
  Function source(mesh, f);
//  FacetNormal n(mesh);
  AvgMeshSize h(mesh);

  // Define PDE
  PoissonBilinearForm a(mesh);
  PoissonLinearForm L(source);
  LinearPDE pde(a, L, mesh);

  // Solve PDE
  Function u(mesh);
  pde.solve(u);

  // Project solution onto continuous basis for post-processing
  Function u_proj(mesh);
  P1ProjectionBilinearForm a_proj(mesh);
  P1ProjectionLinearForm L_proj(u);
  LinearPDE pde_proj(a_proj, L_proj, mesh);
  pde_proj.solve(u_proj);

  // Save solution to file
  File file("poisson.pvd");
  file << u_proj;

  return 0;
}
