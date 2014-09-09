// Copyright (C) 2004-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2006, 2007
//
// First added:  2004
// Last changed: 2007-12-28

#include <dolfin/la/Matrix.h>
#include <dolfin/fem/BoundaryCondition.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/DiscreteFunction.h>
#include <dolfin/io/dolfin_io.h>
#include <dolfin/la/LUSolver.h>
#include <dolfin/la/KrylovSolver.h>
#include <dolfin/pde/LinearPDE.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
LinearPDE::LinearPDE(BilinearForm& a, LinearForm& L, Mesh& mesh) :
    a(a),
    L(L),
    mesh(mesh),
    assembler(mesh),
    not_assembled(true)
{
  message("Creating linear PDE.");
}
//-----------------------------------------------------------------------------
LinearPDE::LinearPDE(BilinearForm& a, LinearForm& L, Mesh& mesh,
                     BoundaryCondition& bc) :
    a(a),
    L(L),
    mesh(mesh),
    assembler(mesh),
    not_assembled(true)
{
  message("Creating linear PDE with one boundary condition.");
  bcs.push_back(&bc);
}
//-----------------------------------------------------------------------------
LinearPDE::LinearPDE(BilinearForm& a, LinearForm& L, Mesh& mesh,
                     Array<BoundaryCondition*>& bcs) :
    a(a),
    L(L),
    mesh(mesh),
    assembler(mesh),
    not_assembled(true)
{
  message("Creating linear PDE with %d boundary condition(s).", bcs.size());
  for (uint i = 0; i < bcs.size(); ++i)
  {
    this->bcs.push_back(bcs[i]);
  }
}
//-----------------------------------------------------------------------------
LinearPDE::~LinearPDE()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void LinearPDE::solve(Function& u)
{
  begin("Solving linear PDE.");

  // Create matrix and vector for assembly
  Matrix A;
  Vector b;
  Vector* x = new Vector();
  u.init(mesh, *x, a, 1);

  // Assemble linear system
  assembler.assemble(A, a, not_assembled);
  assembler.assemble(b, L, not_assembled);
  not_assembled = false;

  // Apply boundary conditions
  for (uint i = 0; i < bcs.size(); ++i)
  {
    bcs[i]->apply(A, b, a);
  }

  // Solve linear system
  const std::string solver_type = get("PDE linear solver");
  if (solver_type == "direct")
  {
    if (u.mesh().is_distributed())
	error("Direct Solvers don't work in parallel please use iterative solvers in parallel with dolfin_set(\"PDE linear solver\", \"iterative\"); or run in serial");
    LUSolver solver;
    solver.set("parent", *this);
    solver.solve(A, *x, b);
  }
  else if (solver_type == "iterative")
  {
    KrylovSolver solver(gmres);
    solver.set("parent", *this);
    solver.solve(A, *x, b);
  }
  else
  {
    error("Unknown solver type \"%s\".", solver_type.c_str());
  }

  end();
}
//-----------------------------------------------------------------------------
void LinearPDE::solve(Function& u0, Function& u1)
{
  // Solve system
  Function u(mesh);
  solve(u);

  // Extract sub functions
  u0 = u[0];
  u1 = u[1];
}
//-----------------------------------------------------------------------------
void LinearPDE::solve(Function& u0, Function& u1, Function& u2)
{
  // Solve system
  Function u(mesh);
  solve(u);

  // Extract sub functions
  u0 = u[0];
  u1 = u[1];
  u2 = u[2];
}
//-----------------------------------------------------------------------------
