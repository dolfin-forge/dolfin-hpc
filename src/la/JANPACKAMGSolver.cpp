// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_JANPACK

#ifdef HAVE_MPI
#include <dolfin/main/MPI.h>
#endif

#include <janpack/amg_solver.h>

#include <dolfin/la/JANPACKMat.h>
#include <dolfin/la/JANPACKVec.h>
#include <dolfin/la/JANPACKAMGSolver.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
JANPACKAMGSolver::JANPACKAMGSolver()
{
}
//-----------------------------------------------------------------------------
JANPACKAMGSolver::~JANPACKAMGSolver()
{
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKAMGSolver::solve(const JANPACKMat& A, JANPACKVec& x, const JANPACKVec& b)
{
  // Check dimensions
  uint M = A.size(0);
  uint N = A.size(1);
  if ( N != b.size() )
    error("Non-matching dimensions for linear system.");
  
  // Write a message
  message("Solving linear system of size %d x %d (AMG solver).", M, N);

  // Reinitialize solution vector if necessary
  if (x.local_size() != b.local_size())
    x.init(b.local_size());


  int num_iterations;
  num_iterations = jp_amg_solver(A.mat(), x.vec(), b.vec(), 
				 JP_AMG_VCYCLE, 0.2, //get("AMG theta"),
				 75,100,1e-3);
  //				 get("AMG maximum iterations"),
  //				 get("AMG relative tolerance"));

  message("AMG solver converged in %d iterations.", num_iterations);
  return num_iterations;
}
//-----------------------------------------------------------------------------
#endif
