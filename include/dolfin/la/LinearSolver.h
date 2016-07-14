// Copyright (C) 2004-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2006.
// Modified by Ola Skavhaug 2008.
//
// First added:  2004-06-19
// Last changed: 2008-05-10

#ifndef __DOLFIN_LINEAR_SOLVER_H
#define __DOLFIN_LINEAR_SOLVER_H

#include <dolfin/parameter/Parametrized.h>

#include <dolfin/la/SolverType.h>
#include <dolfin/la/PreconditionerType.h>

namespace dolfin
{

class GenericMatrix;
class GenericVector;
class KrylovSolver;
class LUSolver;

/// This class provides a general solver for linear systems Ax = b.
/// Available methods are defined in SolverType.h and available
/// preconditioners are defined in PreconditionerType.h.

class LinearSolver : public Parametrized
{

public:

  /// Create linear solver
  LinearSolver(SolverType solver_type = lu, PreconditionerType pc_type = ilu);

  /// Destructor
  ~LinearSolver();

  /// Solve linear system Ax = b
  uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b);

private:

  // LU solver
  LUSolver * lu_solver;

  // Krylov solver
  KrylovSolver * krylov_solver;

};

} /* namespace dolfin */

#endif /* __DOLFIN_LINEAR_SOLVER_H */
