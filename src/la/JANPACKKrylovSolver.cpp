// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_JANPACK

#include <janpack/krylov_solver.h>

#include <dolfin/la/JANPACKMat.h>
#include <dolfin/la/JANPACKVec.h>
#include <dolfin/la/JANPACKKrylovSolver.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
JANPACKKrylovSolver::JANPACKKrylovSolver(SolverType method, 
					 PreconditionerType pc) :
  method(method), pc_janpack(pc)
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

  jp_pc_t pc_type;
  if (pc_janpack == none)
    pc_type = JP_PC_NONE;
  else if(pc_janpack == jacobi)
    pc_type = JP_PC_JACOBI;
  else
    pc_type = JP_PC_NONE;

  int num_iterations;
  num_iterations = jp_krylov_solver(A.mat(), x.vec(), b.vec(), 
				    (jp_solver_t) getType(method), pc_type,
				    get("Krylov maximum iterations"),
				    get("Krylov relative tolerance"));
  message("Krylov solver converged in %d iterations.", num_iterations);
  return num_iterations;
}
//-----------------------------------------------------------------------------
int JANPACKKrylovSolver::getType(SolverType method) const
{

  switch (method)
  {
  case bicgstab:
    return JP_BICGSTAB;
  case cg:
    return JP_CG;
  default:
    warning("Requested Krylov method unknown. Using BICGSTAB.");
    return JP_BICGSTAB;
  }
}
//-----------------------------------------------------------------------------
#endif
