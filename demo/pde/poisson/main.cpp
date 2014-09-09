// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-02-07
// Last changed: 2007-08-20
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
//     u(x, y)     = 0               for x = 0
//     du/dn(x, y) = 25 sin(5 pi y)  for x = 1
//     du/dn(x, y) = 0               otherwise

#include <dolfin/config/dolfin_config.h>
#ifdef ENABLE_UFL 
#include "ufc2/Poisson.h"
#else
#include "ufc1/Poisson.h"
#endif

#include <dolfin/common/constants.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/DirichletBC.h>

#include <dolfin/function/Function.h>
#include <dolfin/mesh/UnitSquare.h>
#include <dolfin/pde/LinearPDE.h>
#include <dolfin/parameter/parameters.h>

using namespace dolfin;

int main()
{
  // Source term
  class Source : public Function
  {
  public:
    
    Source(Mesh& mesh) : Function(mesh) {}
    
    void eval(real * value, const real* x) const
    {
      real dx = x[0] - 0.5;
      real dy = x[1] - 0.5;
      value[0] = 500.0*exp(-(dx*dx + dy*dy)/0.02);
    }

    uint rank() const
    {
      return 0;
    }

    uint dim(uint i) const
    {
      return 1;
  }

  };

  // Neumann boundary condition
  class Flux : public Function
  {
  public:

    Flux(Mesh& mesh) : Function(mesh) {}

    void eval(real * value, const real* x) const
    {
      if (x[0] > DOLFIN_EPS)
        value[0] = 25.0*sin(5.0*DOLFIN_PI*x[1]);
      else
        value[0] = 0.0;
    }

    uint rank() const
    {
      return 0;
    }

    uint dim(uint i) const
    {
      return 1;
    }

  };

  // Sub domain for Dirichlet boundary condition
  class DirichletBoundary : public SubDomain
  {
    bool inside(const real* x, bool on_boundary) const
    {
      return x[0] < DOLFIN_EPS && on_boundary;
    }
  };

  // Create mesh
  //UnitSquare mesh(32, 32);
  Mesh mesh("UnitSquareMesh_32x32.xml");

  // Create functions
  Source f(mesh);
  Flux g(mesh);

  // Create boundary condition
  Function u0(mesh, 0.0);
  DirichletBoundary boundary;
  DirichletBC bc(u0, mesh, boundary);
  
  // Define PDE
  PoissonBilinearForm a(mesh);
  PoissonLinearForm L(f, g);
  dolfin_set("PDE linear solver", "iterative");
  LinearPDE pde(a, L, mesh, bc);

  // Solve PDE
  Function u(mesh);

  pde.solve(u);
   
  // Save solution to file
  File file("poisson.pvd");
  file << u;

  return 0;
}
