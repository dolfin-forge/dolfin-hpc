// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// This demo solves the Stokes equations, using quadratic
// elements for the velocity and first degree elements for
// the pressure (Taylor-Hood elements). The sub domains
// for the different boundary conditions used in this
// simulation are computed by the demo program in
// src/demo/mesh/subdomains.

#include "Stokes.h"
#include "../Stokes2D.h"

#include <dolfin.h>

using namespace dolfin;

int main( int argc, char ** argv )
{
	dolfin_init( argc, argv );

	{
		// Read mesh and sub domain markers
		Mesh mesh( "../../../../data/meshes/dolfin-2.bin" );

		// Create functions for boundary conditions
		Analytic< NoSlip >    noslip( mesh );
		Analytic< Inflow >    inflow( mesh );
		Analytic< Zero< 1 > > zero( mesh );

		// Define sub systems for boundary conditions
		SubSystem velocity( 0 );
		SubSystem pressure( 1 );

		// No-slip boundary condition for velocity
		NoSlipBoundary nsb;
		DirichletBC    bc0( noslip, mesh, nsb, velocity );

		// Inflow boundary condition for velocity
		RightBoundary rb;
		DirichletBC   bc1( inflow, mesh, rb, velocity );

		// Boundary condition for pressure at outflow
		LeftBoundary lb;
		DirichletBC  bc2( zero, mesh, lb, pressure );

		// Set up PDE
		Analytic< Zero< 2 > >    f( mesh );
		Stokes::BilinearForm a( mesh );
		Stokes::LinearForm   L( mesh, f );

		// Solve PDE
		Matrix A;
		Vector b;

		a.assemble( A, true );
		L.assemble( b, true );
		bc0.apply( A, b, a );
		bc1.apply( A, b, a );
		bc2.apply( A, b, a );

		Function     u( a.trial_space() );
		Function     p( a.trial_space() );
		KrylovSolver solver( bicgstab, bjacobi );

		solver.solve( A, u.vector(), b );
		solver.solve( A, p.vector(), b );
		u.sync();

		// Save solution in VTK format
		File ufile_pvd( "velocity.pvd" );
		ufile_pvd << u;
		File pfile_pvd( "pressure.pvd" );
		pfile_pvd << p;
	}

	dolfin_finalize();

	return 0;
}
