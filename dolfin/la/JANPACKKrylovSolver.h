// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __JANPACK_KRYLOV_SOLVER_H
#define __JANPACK_KRYLOV_SOLVER_H

#ifdef HAS_JANPACK
#include <dolfin/common/types.h>
#include <dolfin/parameter/Parametrized.h>
#include "SolverType.h"
#include "PreconditionerType.h"

namespace dolfin 
{
  /// Forward declarations
  class JANPACKMat;
  class JANPACKVec;

  class JANPACKKrylovSolver : public Parametrized
  {
  public:

    /// Create Krylov solver for a particular method and preconditioner
    JANPACKKrylovSolver(SolverType method=default_solver, PreconditionerType pc=default_pc);

    /// Solve linear system Ax = b and return number of iterations
    uint solve(const JANPACKMat& A, JANPACKVec& x, const JANPACKVec& b); 
    
  private:

    /// Krylov method
    SolverType method;

  };
}

#endif

#endif
