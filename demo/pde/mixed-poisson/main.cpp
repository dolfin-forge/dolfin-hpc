// Copyright (C) 2007 Anders Logg and Marie Rognes.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2008.
//
// First added:  2007-04-20
// Last changed: 2008-07-12
//
// This demo program solves the mixed formulation of
// Poisson's equation:
//
//     sigma + grad(u) = 0
//          div(sigma) = f
//
// The corresponding weak (variational problem)
//
//     <tau, sigma> - <div(tau), u> = 0       for all tau
//                  <w, div(sigma)> = <w, f>  for all w
//
// is solved using BDM (Brezzi-Douglas-Marini) elements
// of degree q (tau, sigma) and DG (discontinuous Galerkin)
// elements of degree q - 1 for (w, u).


#include "ufc1/MixedPoisson.h"
#include "ufc1/P1Projection.h"

#include <dolfin/common/constants.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/function/Function.h>
#include <dolfin/mesh/UnitSquare.h>
#include <dolfin/pde/LinearPDE.h>

using namespace dolfin;

int main()
{
  // Source term
  class Source : public Function
  {
  public:
    
    Source(Mesh& mesh) : Function(mesh) {}

    void eval(real * value, const real* x) const
    {
      real dx = x[0] - 0.5;
      real dy = x[1] - 0.5;
      value[0] = 500.0*exp(-(dx*dx + dy*dy)/0.02);
    }

    uint rank() const
    {
      return 0;
    }

    uint dim(uint i) const
    {
      return 1;
    }

  };

  // Create mesh and source term
  UnitSquare mesh(16, 16);
  Source f(mesh);
  
  // Define PDE
  MixedPoissonBilinearForm a;
  MixedPoissonLinearForm L(f);
  LinearPDE pde(a, L, mesh);

  // Solve PDE
  Function sigma;
  Function u;
  pde.solve(sigma, u);

  // Project sigma onto P1 continuous Lagrange for post-processing
  Function sigma_projected;
  P1ProjectionBilinearForm a_projection;
  P1ProjectionLinearForm L_projection(sigma);
  LinearPDE pde_project(a_projection, L_projection, mesh);
  pde_project.solve(sigma_projected);


  // Save solution to file
  File f0("sigma.xml");
  File f1("u.xml");
  f0 << sigma;
  f1 << u;

  // Save solution to pvd format
  File f3("sigma.pvd");
  File f4("u.pvd");
  f3 << sigma_projected;
  f4 << u;

  return 0;
}
