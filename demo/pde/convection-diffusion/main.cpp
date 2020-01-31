// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// This demo solves the time-dependent convection-diffusion equation by
// a least-squares stabilized cG(1)cG(1) method. The velocity field used
// in the simulation is the output from the Stokes (Taylor-Hood) demo.
// The sub domains for the different boundary conditions are computed
// by the demo program in src/demo/subdomains.

#include "ConvectionDiffusion.h"
#include <dolfin.h>

using namespace dolfin;

// Sub domain for Dirichlet boundary condition
struct DirichletBoundary : public SubDomain
{
	bool inside( const real * x, bool on_boundary ) const
	{
		std::cout << x[0] << " " << x[1] << std::endl;
		return x[0] > .99;
	}

	// MeshValues< uint, Cell > sub_domains;
};

int main( int argc, char * argv[] )
{
	dolfin_init( argc, argv );
	// Read mesh and sub domain markers
	Mesh                     mesh( "../../../data/meshes/dolfin-2.xml.gz" );
	MeshValues< uint, Cell > sub_domains( mesh );
	File( "subdomains.xml.gz" ) >> sub_domains;

	// Convection velocity, source term and initial condition
	Function velocity( mesh );
	Function f( mesh );
	Function u0( mesh );

	// Set up forms
	ConvectionDiffusion::BilinearForm a( mesh, velocity );
	ConvectionDiffusion::LinearForm   L( mesh, u0, velocity, f );
	velocity.init(L.test_space() );
	File( "velocity.bin" ) >> velocity;

	// Create finite element spaces
	FiniteElementSpace * FE_vel  = L.create_coefficient_space( "b" );
	FiniteElementSpace * FE_scal = L.create_coefficient_space( "u0" );

	// Read the velocity from file
	velocity.init( *FE_vel );

	// Set up boundary condition
	Function          g( a.trial_space() );
	Function          g0( a.trial_space() );
	DirichletBoundary boundary;
	DirichletBC       bc( g, mesh, boundary );
	DirichletBC       bc0( g0, mesh, boundary );

	// Linear system
	Matrix A;
	Vector b;

	// Solution vector
	Function u1( *FE_scal );

	// LU
	LUSolver lu;

	// Assemble matrix
	Assembler::assemble( A, a, true );
	Assembler::assemble( b, L, true );

	real T = 2.0;
	real k = 0.05;
	real t = k;

	// Output file
	File file( "temperature.pvd" );

	// Time-stepping
	while ( t < T )
	{
		// Assemble vector and apply boundary conditions
		Assembler::assemble( b, L, false );
		bc.apply( A, b, a );
		bc0.apply( A, b, a );

		// Solve the linear system
		// lu.factorized_solve(u1.vector(), b);
		lu.solve( A, u1.vector(), b );

		// Save the solution to file
		file << u1;

		// Move to next interval
		t += k;
		u0 = u1;
	}

	delete FE_vel;
	delete FE_scal;

	dolfin_finalize();
	return 0;
}
