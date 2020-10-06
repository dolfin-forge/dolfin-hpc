// Copyright (C) 2006-2007 Anders Logg and Kristian Oelgaard.
// Licensed under the GNU LGPL Version 2.1.
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

#include "P1Projection.h"
#include "Poisson.h"

using namespace dolfin;

// Source term
struct Source : public Value< Source, 1 >
{
  void eval( real * values, const real * x ) const
  {
    real dx   = x[0] - 0.5;
    real dy   = x[1] - 0.5;
    values[0] = 500.0 * exp( -( dx * dx + dy * dy ) / 0.02 );
  }
};

int main()
{

  dolfin_init();
  
  // Create mesh
  UnitSquare mesh( 64, 64 );

  // Create functions
  Analytic< Source > f( mesh );

  // Define PDE
  Poisson::BilinearForm a( mesh );
  Poisson::LinearForm   L( mesh, f );

  // Solve PDE
  Function u( a.trial_space() );
  Matrix   A;
  Vector   b;
  a.assemble( A, true );
  L.assemble( b, true );

  KrylovSolver solver( bicgstab, bjacobi );

  solver.solve( A, u.vector(), b );
  u.sync();

  // Project solution onto continuous basis for post-processing
  P1Projection::BilinearForm a_proj( mesh );
  P1Projection::LinearForm   L_proj( mesh, u );
  Function                   u_proj( a_proj.trial_space() );

  Matrix A_proj;
  Vector b_proj;
  a_proj.assemble( A_proj, true );
  L_proj.assemble( b_proj, true );
  solver.solve( A_proj, u_proj.vector(), b_proj );
  u_proj.sync();

  // Save solution to file
  File file( "poisson.pvd" );
  file << u_proj;

  dolfin_finalize();
  
  return 0;
}
