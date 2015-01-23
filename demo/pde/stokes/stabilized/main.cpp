// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-02-09
// Last changed: 2007-07-11
//
// This demo solves the Stokes equations, using stabilized
// first order elements for the velocity and pressure. The
// sub domains for the different boundary conditions used
// in this simulation are computed by the demo program in
// src/demo/mesh/subdomains.

#include <dolfin.h>

#include "../Stokes2D.h"
#include "Stokes.h"

using namespace dolfin;

int main(int argc, char* argv[])
{
  dolfin_init(argc, argv);

	{
		// Read mesh and sub domain markers
		Mesh mesh("../../../../data/meshes/dolfin-2.xml.gz");

		// Create functions for boundary conditions
		NoSlip nslip;
		Function noslip(mesh, nslip);
		Inflow in;
		Function inflow(mesh, in);
		Function zero(mesh, 0.0);

		// Define sub systems for boundary conditions
		SubSystem velocity(0);
		SubSystem pressure(1);

		// No-slip boundary condition for velocity
		NoSlipBoundary nsb;
		DirichletBC bc0(noslip, mesh, nsb, velocity);

		// Inflow boundary condition for velocity
		RightBoundary rb;
		DirichletBC bc1(inflow, mesh, rb, velocity);

		// Boundary condition for pressure at outflow
		LeftBoundary lb;
		DirichletBC bc2(zero, mesh, lb, pressure);

		// Collect boundary conditions
		Array <BoundaryCondition*> bcs;
		bcs.push_back(&bc0);
  	bcs.push_back(&bc1);
	  bcs.push_back(&bc2);

		// Set up PDE
		MeshSize h(mesh);
		Function f(mesh, 2, 0.0);
		StokesBilinearForm a(h);
		StokesLinearForm L(f, h);
		LinearPDE pde(a, L, mesh, bcs);

		// Solve PDE
		Function u(mesh);
		Function p(mesh);
		pde.set("PDE linear solver", "direct");
		pde.solve(u, p);

		// Save solution
		File ufile("velocity.xml");
		ufile << u;
		File pfile("pressure.xml");
		pfile << p;

		// Save solution in VTK format
		File ufile_pvd("velocity.pvd");
		ufile_pvd << u;
		File pfile_pvd("pressure.pvd");
		pfile_pvd << p;

		File x_file("x.xml");
		x_file << u.vector();
	}

	dolfin_finalize();
}
