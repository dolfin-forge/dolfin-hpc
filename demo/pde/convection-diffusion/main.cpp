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
#ifdef ENABLE_UFL 
#include "ufc2/ConvectionDiffusion.h"
#else
#include "ufc1/ConvectionDiffusion.h"
#endif

using namespace dolfin;

int main(int argc, char *argv[])
{
  // Read mesh and sub domain markers
  Mesh mesh("../../../data/meshes/dolfin-2.xml.gz");
  MeshFunction<unsigned int> sub_domains(mesh, "subdomains.xml.gz");

  // Read velocity field
  Function velocity(mesh);
#ifdef ENABLE_UFL 
  File vel("ufc2/velocity.xml.gz");
  vel >> velocity;
#else
  File vel("ufc1/velocity.xml.gz");
  vel >> velocity;
#endif

  // Source term and initial condition
  Function f(mesh, 0.0);
  Function u0(mesh, 0.0);

  // Set up forms
  ConvectionDiffusionBilinearForm a(velocity);
  ConvectionDiffusionLinearForm L(u0, velocity, f);

  // Set up boundary condition
  Function g(mesh, 1.0);
  DirichletBC bc(g, sub_domains, 1);

  // Linear system
  Matrix A;
  Vector x, b;

  // Solution vector
  Function u1(mesh, a, 1);

  // LU
  LUSolver lu;

  // Assemble matrix
  Assembler assembler(mesh);
  assembler.assemble(A, a, true);
  assembler.assemble(b, L, true);
  bc.apply(A, b, a);
  //lu.factorize(A);

  // Parameters for time-stepping
  real T = 2.0;
  real k = 0.05;
  real t = k;
  
  // Output file
  File file("temperature.pvd");

  // Time-stepping
#ifndef NO_PROGRESS_BAR
  Progress p("Time-stepping");
#endif
  while ( t < T )
  {
    // Assemble vector and apply boundary conditions
    assembler.assemble(b, L, false);
    bc.apply(A, b, a);
    
    // Solve the linear system
    //lu.factorized_solve(x, b);
    lu.solve(A, x, b);
    
    // Save the solution to file
    file << u1;

    // Move to next interval
#ifndef NO_PROGRESS_BAR
    p = t / T;
#endif
    t += k;
    u0 = u1;
  }


}
