// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#ifdef HAS_JANPACK

#include <cg.h>

#include "JANPACKMat.h"
#include "JANPACKVec.h"
#include "JANPACKKrylovSolver.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
JANPACKKrylovSolver::JANPACKKrylovSolver(SolverType method, PreconditionerType pc)
{
}
//-----------------------------------------------------------------------------
uint JANPACKKrylovSolver::solve(const JANPACKMat& A, JANPACKVec& x, const JANPACKVec& b)
{
  // Check dimensions
  uint M = A.size(0);
  uint N = A.size(1);
  if ( N != b.size() )
    error("Non-matching dimensions for linear system.");
  
  // Write a message
  message("Solving linear system of size %d x %d (Krylov solver).", M, N);

  // Reinitialize solution vector if necessary
  x.init(M);

  int num_iterations = cg_crs(A.mat(), x.vec(), b.vec());
  message("Krylov solver converged in %d iterations.", num_iterations);
  return 50;
}
//-----------------------------------------------------------------------------

#endif
