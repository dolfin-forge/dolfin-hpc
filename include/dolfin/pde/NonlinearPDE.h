// Copyright (C) 2005-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2007
//
// First added:  2005-10-24
// Last changed: 2007-05-15

#ifndef __NONLINEAR_PDE_H
#define __NONLINEAR_PDE_H

#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/BilinearForm.h>
#include <dolfin/fem/LinearForm.h>
#include <dolfin/nls/NonlinearProblem.h>
#include <dolfin/nls/NewtonSolver.h>

namespace dolfin
{

/// This class provides automated solution of nonlinear PDEs.

class NonlinearPDE : public NonlinearProblem, public Parametrized
{
public:

  /// Constructor
  NonlinearPDE(BilinearForm& a, LinearForm& L, Mesh& mesh,
               BoundaryCondition& bc);

  /// Constructor
  NonlinearPDE(BilinearForm& a, LinearForm& L, Mesh& mesh,
               Array<BoundaryCondition*>& bcs);

  /// Destructor
  ~NonlinearPDE();

  /// Function called before Jacobian matrix and RHS vector are formed. Users
  /// can supply this function to perform updates.
  virtual void update(GenericVector const& x);

  /// User-defined function to compute F(u) its Jacobian
  void form(GenericMatrix& A, GenericVector& b, GenericVector const& x);

  /// Solve PDE
  void solve(Function& u, real& t, real const& T, real const& dt);

private:

  // The bilinear form
  BilinearForm& a;

  // The linear form
  LinearForm& L;

  // The mesh
  Mesh& mesh;

  // The boundary conditions
  Array<BoundaryCondition*> bcs;

  // The solution vector
  Vector x;

  // Assembler
  Assembler assembler;

  // Solver
  NewtonSolver newton_solver;

  bool reset;
};

}

#endif
