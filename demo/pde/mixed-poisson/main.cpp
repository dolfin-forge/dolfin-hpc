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

// DOES NOT WORK YET!!!
// FIXME: add support for ufl mixed elements

#include "MixedPoisson.h"
#include "P1Projection.h"

#include <dolfin.h>

using namespace dolfin;

// Source term
struct Source : public Value<Source,1>
{
  void eval( real * value, const real * x ) const
  {
    real dx  = x[0] - 0.5;
    real dy  = x[1] - 0.5;
    value[0] = 500.0 * exp( -( dx * dx + dy * dy ) / 0.02 );
  }
};

int main()
{
	// Create mesh and source term
	UnitSquare mesh( 16, 16 );
  Analytic<Source>  f(mesh);

	// Solve PDE
    MixedPoisson::BilinearForm a( mesh );
    MixedPoisson::LinearForm   L( mesh, f );

    // Solve PDE
    Matrix A;
    Vector b;

    a.assemble(A, true);
    L.assemble(b, true);

    Function sigma(a.trial_space());
    KrylovSolver solver(bicgstab, bjacobi);

    solver.solve(A, sigma.vector(), b);
    sigma.sync();

    // Function u( u[0] );
    // Function sigma( u[1] );

    // Save solution to file
    File( "sigma.pvd" ) << sigma;

	// Project sigma onto P1 continuous Lagrange for post-processing
    P1Projection::BilinearForm a_projection( mesh );
    P1Projection::LinearForm   L_projection( mesh, sigma );

    // Solve PDE
    Matrix A_proj;
    Vector b_proj;

    a_projection.assemble(A_proj, true);
    L_projection.assemble(b_proj, true);

    Function sigma_projected(a_projection.trial_space());
    KrylovSolver solver_proj(bicgstab, bjacobi);

    solver_proj.solve(A_proj, sigma_projected.vector(), b_proj);
    sigma_projected.sync();

    // Save solution to file
    File( "sigma_projected.pvd" ) << sigma_projected;

	return 0;
}
