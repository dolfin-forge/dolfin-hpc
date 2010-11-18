// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_JANPACK

#include <janpack/bicgstab.h>
#include <janpack/cg.h>

#include "JANPACKMat.h"
#include "JANPACKVec.h"
#include "JANPACKKrylovSolver.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
JANPACKKrylovSolver::JANPACKKrylovSolver(SolverType method, 
					 PreconditionerType pc) :
  method(method)
{
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKKrylovSolver::solve(const JANPACKMat& A, JANPACKVec& x, const JANPACKVec& b)
{
  // Check dimensions
  uint M = A.size(0);
  uint N = A.size(1);
  if ( N != b.size() )
    error("Non-matching dimensions for linear system.");
  
  // Write a message
  message("Solving linear system of size %d x %d (Krylov solver).", M, N);

  // Reinitialize solution vector if necessary
  x.init(b.local_size());

  int num_iterations;
  switch (method)
  {
  case bicgstab:
    num_iterations = bicgstab_crs(A.mat(), x.vec(), b.vec()); break;
  case cg:
    num_iterations = cg_crs(A.mat(), x.vec(), b.vec()); break;
  default:
    error("Krylov solver not supported by JANPACK");
  }  
    
  message("Krylov solver converged in %d iterations.", num_iterations);
  return num_iterations;
}
//-----------------------------------------------------------------------------

#endif
