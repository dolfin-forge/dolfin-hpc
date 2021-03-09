// Copyright (C) 2007 Kristian B. Oelgaard.
// Licensed under the GNU LGPL Version 2.1.
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
struct DirichletBoundary : public SubDomain
{
  bool inside( const real * x, bool on_boundary ) const
  {
    // return ( std::abs( x[0] ) < DOLFIN_EPS );
    return x[0] < DOLFIN_EPS && on_boundary;
  }
};

// Source term
struct Source : public Value< Source, 1 >
{
  void eval( real * values, const real * x ) const
  {
    values[0] = 9.0 * DOLFIN_PI * DOLFIN_PI * sin( 3.0 * DOLFIN_PI * x[0] );
  }
};

int main()
{
  dolfin_init();

  // Create mesh
  UnitInterval mesh( 50 );

  // Set up BCs
  Constant          zero( 0.0 );
  DirichletBoundary boundary;
  DirichletBC       bc( zero, mesh, boundary );

  // Create source
  Analytic< Source > source( mesh );

  // Define PDE
  Poisson::BilinearForm a( mesh );
  Poisson::LinearForm   L( mesh, source );

  // Solve PDE
  Matrix A;
  Vector b;

  a.assemble( A, true );
  L.assemble( b, true );
  bc.apply( A, b, a );

  Function u( a.trial_space() );

  KrylovSolver solver( bicgstab, bjacobi );
  solver.solve( A, u.vector(), b );

  // Save solution to file
  File( "poisson.pvd" ) << u;

  dolfin_finalize();

  return 0;
}
