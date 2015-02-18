// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-02-09
// Last changed: 2007-07-11
//
// This demo solves the time-dependent convection-diffusion equation by
// a least-squares stabilized cG(1)cG(1) method. The velocity field used
// in the simulation is the output from the Stokes (Taylor-Hood) demo.
// The sub domains for the different boundary conditions are computed
// by the demo program in src/demo/subdomains.

#include <dolfin.h>
#include "ConvectionDiffusion.h"

using namespace dolfin;

int main(int argc, char *argv[]) 
{
  dolfin_init(argc, argv);
  // Read mesh and sub domain markers
  Mesh mesh("../../../data/meshes/dolfin-2.xml.gz");
  MeshFunction<unsigned int> sub_domains(mesh, "subdomains.xml.gz");

  // Convection velocity, source term and initial condition
  Function velocity(mesh);
  Function f(mesh, 0.0);
  Function u0(mesh, 0.0);

  // Set up forms
  ConvectionDiffusionBilinearForm a(velocity);
  ConvectionDiffusionLinearForm L(u0, velocity, f);

  // Create finite element spaces
  FiniteElementSpace * FE_vel = L.create_coefficient_space("b");
  FiniteElementSpace * FE_scal = L.create_coefficient_space("u0");

  // Read the velocity from file
  velocity.init(*FE_vel);
  File vel("velocity.xml.gz");
  vel >> velocity;

  // Set up boundary condition
  Function g(mesh, 1.0);
  Function g0(mesh, 0.0);
  DirichletBC bc(g, sub_domains, 1);
  DirichletBC bc0(g0, sub_domains, 0);

  // Linear system
  Matrix A;
  Vector b;

  // Solution vector
  Function u1(*FE_scal);

  // LU
  LUSolver lu;

  // Assemble matrix
  Assembler assembler(mesh);
  assembler.assemble(A, a, true);
  assembler.assemble(b, L, true);

  real T = 2.0;
  real k = 0.05;
  real t = k;

  // Output file
  File file("temperature.pvd");

  // Time-stepping
  Progress p("Time-stepping");
  while (t < T) 
    {
      // Assemble vector and apply boundary conditions
      assembler.assemble(b, L, false);
      bc.apply(A, b, a);
      bc0.apply(A, b, a);

      // Solve the linear system
      //lu.factorized_solve(u1.vector(), b);
      lu.solve(A, u1.vector(), b);

      // Save the solution to file
      file << u1;

      // Move to next interval
      p = t / T;
      t += k;
      u0 = u1;
    }

  delete FE_vel, FE_scal;

  dolfin_finalize();
  return 0;
}
