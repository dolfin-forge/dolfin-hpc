// Copyright (C) 2006-2007 Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2008
//
// First added:  2006-02-07
// Last changed: 2008-05-21
//
// This demo program solves the equations of static
// linear elasticity for a gear clamped at two of its
// ends and twisted 30 degrees.

#include "Elasticity.h"
#include <dolfin.h>
#include <dolfin/config/dolfin_config.h>

using namespace dolfin;

struct ConstantFunction : public Value<ConstantFunction,2> //Value< ConstantFunction, 3 >
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

struct ConstantFunction1 : public Value<ConstantFunction1,1,3> //Value< ConstantFunction, 3 >
{
  ConstantFunction1( real c_ )
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
struct Clamp : public Value< Clamp >
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
struct Rotation : public Value< Rotation >
{
	void eval( real * values, const real * x ) const
	{
		// Center of rotation
		real y0 = 0.5;
		real z0 = 0.219;

		// Angle of rotation (30 degrees)
		real theta = 0.5236;

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
	// Function f( mesh, 3, 0.0 );
	ConstantFunction             f_( 0.0 );
	Analytic< ConstantFunction > f( mesh, f_ );

	// Set up boundary condition at left end
	// Clamp       c;
	// Function    clamp( mesh, c );
	Analytic< Clamp > clamp( mesh );
	Left              left;
	DirichletBC       bcl( clamp, mesh, left );

	// Set up boundary condition at right end
	// Rotation    r;
	// Function    rotation( mesh, r );
	Analytic< Rotation > rotation( mesh );
	Right                right;
	DirichletBC          bcr( rotation, mesh, right );

	// Set up boundary conditions
	// Array< BoundaryCondition * > bcs;
	// bcs.push_back( &bcl );
	// bcs.push_back( &bcr );

	// Set elasticity parameters
	real                         E  = 10.0;
	real                         nu = 0.3;
	ConstantFunction1             mu_( E / ( 2 * ( 1 + nu ) ) );
	Analytic< ConstantFunction1 > mu( mesh, mu_ );
	ConstantFunction1 lambda_( E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) ) );
	Analytic< ConstantFunction1 > lambda( mesh, lambda_ );

	// Set up PDE
	Elasticity::BilinearForm   a( mesh, mu, lambda );
	Elasticity::LinearForm     L( mesh, f );

	// Solve PDE
	Matrix A;
	Vector b;

	a.assemble( A, true );
	L.assemble( b, true );
	bcl.apply( A, b, a );
	bcr.apply( A, b, a );

	Function     u( a.trial_space() );
	KrylovSolver solver( bicgstab, bjacobi );

	solver.solve( A, u.vector(), b );
	u.sync();

	// Solve PDE (using direct solver)
	// Function u( a.trial_space() );
	// pde.set( "PDE linear solver", "direct" );
	// pde.solve( u );

	// Save solution to VTK format
	File vtk_file( "elasticity.pvd" );
	vtk_file << u;

	// Save solution to XML format
	File xml_file( "elasticity.xml" );
	xml_file << u;

	return 0;
}
