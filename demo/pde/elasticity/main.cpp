// Copyright (C) 2006-2007 Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// This demo program solves the equations of static
// linear elasticity for a gear clamped at two of its
// ends and twisted 30 degrees.

#include "Elasticity.h"

#include <dolfin.h>

using namespace dolfin;

struct Null : public Value< Null, 3 >
{
  void eval( real * values, const real * x ) const
  {
    values[0] = 0.;
    values[1] = 0.;
    values[2] = 0.;
  }
};

struct ConstantFunction : public Value< ConstantFunction, 1, 3 >
{
  ConstantFunction( real c_ )
    : c( c_ )
  {
  }

  void eval( real * values, const real * x ) const
  {
    values[0] = c;
    values[1] = c;
    values[2] = c;
  }

  real c;
};

// Dirichlet boundary condition for clamp at left end
struct Clamp : public Value< Clamp, 3 >
{
  void eval( real * values, const real * x ) const
  {
    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;
  }
};

// Sub domain for clamp at left end
class Left : public SubDomain
{
  bool inside( const real * x, bool on_boundary ) const
  {
    return x[0] < 0.5 && on_boundary;
  }
};

// Dirichlet boundary condition for rotation at right end
struct Rotation : public Value< Rotation, 3 >
{
  void eval( real * values, const real * x ) const
  {
    // Center of rotation
    real const y0 = 0.5;
    real const z0 = 0.219;

    // Angle of rotation (30 degrees)
    real const theta = 0.5236;

    // New coordinates
    real y = y0 + ( x[1] - y0 ) * cos( theta ) - ( x[2] - z0 ) * sin( theta );
    real z = z0 + ( x[1] - y0 ) * sin( theta ) + ( x[2] - z0 ) * cos( theta );

    // Clamp at right end
    values[0] = 0.0;
    values[1] = y - x[1];
    values[2] = z - x[2];
  }
};

// Sub domain for rotation at right end
class Right : public SubDomain
{
  bool inside( const real * x, bool on_boundary ) const
  {
    return x[0] > 0.9 && on_boundary;
  }
};

int main()
{
  // Read mesh
  Mesh mesh( "../../../data/meshes/gear.xml.gz" );

  // Create right-hand side
  Analytic< Null > f( mesh );

  // Set up boundary condition at left end
  Analytic< Clamp > clamp( mesh );
  Left              left;
  DirichletBC       bcl( clamp, mesh, left );

  // Set up boundary condition at right end
  Analytic< Rotation > rotation( mesh );
  Right                right;
  DirichletBC          bcr( rotation, mesh, right );

  // Set elasticity parameters
  real const                   E  = 10.0;
  real const                   nu = 0.3;
  ConstantFunction             mu_( E / ( 2 * ( 1 + nu ) ) );
  Analytic< ConstantFunction > mu( mesh, mu_ );
  ConstantFunction lambda_( E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) ) );
  Analytic< ConstantFunction > lambda( mesh, lambda_ );

  // Set up PDE
  Elasticity::BilinearForm a( mesh, mu, lambda );
  Elasticity::LinearForm   L( mesh, f );

  // Solve PDE
  Matrix A;
  Vector b;

  a.assemble( A, true );
  L.assemble( b, true );
  bcl.apply( A, b, a );
  bcr.apply( A, b, a );

  Function     u( a.trial_space() );
  LinearSolver solver;

  solver.solve( A, u.vector(), b );
  u.sync();

  // Save solution to VTK format
  File vtk_file( "elasticity.pvd" );
  vtk_file << u;

  return 0;
}
