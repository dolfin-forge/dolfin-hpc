// Copyright (C) 2004-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LINEAR_SOLVER_H
#define __DOLFIN_LINEAR_SOLVER_H

#include <dolfin/la/PreconditionerType.h>
#include <dolfin/la/SolverType.h>

namespace dolfin
{

class GenericMatrix;
class GenericVector;
class KrylovSolver;
class LUSolver;

//-----------------------------------------------------------------------------

/// This class provides a general solver for linear systems Ax = b.
/// Available methods are defined in SolverType.h and available
/// preconditioners are defined in PreconditionerType.h.

class LinearSolver
{

public:
  /// Create linear solver
  LinearSolver( SolverType         solver_type = lu,
                PreconditionerType pc_type     = PreconditionerType::ilu );

  /// Destructor
  ~LinearSolver();

  /// Solve linear system Ax = b
  auto solve( const GenericMatrix & A,
              GenericVector &       x,
              const GenericVector & b ) -> size_t;

private:
  // LU solver
  LUSolver * lu_solver { nullptr };

  // Krylov solver
  KrylovSolver * krylov_solver { nullptr };
};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_LINEAR_SOLVER_H */
