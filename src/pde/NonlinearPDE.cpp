// Copyright (C) 2005-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2006-2007.
//
// First added:  2005-10-24
// Last changed: 2007-08-28

#include <dolfin/fem/BoundaryCondition.h>
#include <dolfin/function/Function.h>
#include <dolfin/pde/NonlinearPDE.h>
#include <dolfin/fem/Form.h>
#include <dolfin/log/dolfin_log.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
NonlinearPDE::NonlinearPDE(BilinearForm& a, LinearForm& L, Mesh& mesh,
                           BoundaryCondition& bc) :
    a(a),
    L(L),
    mesh(mesh),
    assembler(mesh),
    reset(true)
{
  message("Creating nonlinear PDE with %d boundary condition(s).", bcs.size());

  // Create array with one boundary condition
  bcs.push_back(&bc);
}
//-----------------------------------------------------------------------------
NonlinearPDE::NonlinearPDE(BilinearForm& a, LinearForm& L, Mesh& mesh,
                           Array<BoundaryCondition*>& bcs) :
    a(a),
    L(L),
    mesh(mesh),
    bcs(bcs),
    assembler(mesh),
    reset(true)
{
  message("Creating nonlinear PDE with %d boundary condition(s).", bcs.size());
}
//-----------------------------------------------------------------------------
NonlinearPDE::~NonlinearPDE()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void NonlinearPDE::update(GenericVector const& x)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void NonlinearPDE::form(GenericMatrix& A, GenericVector& b,
                        GenericVector const& x)
{
  // Assemble
  assembler.assemble(A, a, reset);
  assembler.assemble(b, L, reset);

  reset = false;
  // Apply boundary conditions
  for (uint i = 0; i < bcs.size(); i++)
  {
    bcs[i]->apply(A, b, x, a);
  }
}
//-----------------------------------------------------------------------------
void NonlinearPDE::solve(Function& u, real& t, real const& T, real const& dt)
{
  begin("Solving nonlinear PDE.");

  // Initialise function
  u.init(mesh, x, a, 1);

  // Solve
  while (t < T)
  {
    t += dt;
    newton_solver.solve(*this, x);
  }

  end();
}
//-----------------------------------------------------------------------------
