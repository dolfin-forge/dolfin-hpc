// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/config/dolfin_config.h>

#if defined(HAVE_JANPACK) && !defined(HAVE_JANPACK_MPI)

#include <janpack/amg_solver.h>

#include <dolfin/la/JANPACKMat.h>
#include <dolfin/la/JANPACKVec.h>
#include <dolfin/la/JANPACKAMGSolver.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
JANPACKAMGSolver::JANPACKAMGSolver(MultigridScheme scheme, 
				   MultigridSmoother smoother,
				   MultigridCoarsening cscheme)
  : scheme(scheme), smoother(smoother), cscheme(cscheme), mls_init(false)
{
}
//-----------------------------------------------------------------------------
JANPACKAMGSolver::~JANPACKAMGSolver()
{
  jp_mls_free(mls);
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKAMGSolver::solve(const JANPACKMat& A,
				     JANPACKVec& x, const JANPACKVec& b)
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

  if (mls_init == false) {
    jp_mls_init(mls, b.local_size());
    mls_init = true;
  }

  if (get("AMG keep levels"))
    jp_mls_setopt(mls, JP_ML_KEEP_LEVELS);

  int num_iterations;
  num_iterations = jp_amg_solver(A.mat(), x.vec(), b.vec(), mls,
				 (jp_amg_scheme_t) getScheme(scheme), 
				 (jp_amg_smoother_t) getSmoother(smoother),
				 (jp_amg_cscheme_t) getCoarsening(cscheme),
				 get("AMG theta"), 
				 get("AMG pre-smoothing steps"),
				 get("AMG post-smoothing steps"),
				 get("AMG levels"),
				 get("AMG maximum iterations"),
				 get("AMG relative tolerance"));
  
  message("AMG solver converged in %d iterations.", num_iterations);

  return num_iterations;
}
//-----------------------------------------------------------------------------
int JANPACKAMGSolver::getScheme(MultigridScheme scheme) const
{

  switch (scheme)
  {
  case mgv:
    return JP_AMG_VCYCLE;
  case mgw:
    return JP_AMG_WCYCLE;
  case fmv:
    return JP_AMG_FMV;
  default:
    warning("Requested Multigrid scheme unknown. Using V-cycle.");
    return JP_AMG_VCYCLE;
  }
}
//-----------------------------------------------------------------------------
int JANPACKAMGSolver::getSmoother(MultigridSmoother smoother) const
{

  switch (smoother)
  {
  case mg_jacobi:
    return JP_AMG_JACOBI;
  case mg_gauss_seidel:
    return JP_AMG_GAUSS_SEIDEL;
  case mg_cf_gauss_seidel:
    return JP_AMG_CF_GAUSS_SEIDEL;
  case mg_sor:
    return JP_AMG_SOR;
  case mg_cf_sor:
    return JP_AMG_CF_SOR;
  default:
    warning("Requested Multigrid smoother unknown. Using Gauss-Seidel.");
    return JP_AMG_GAUSS_SEIDEL;
  }
}
//-----------------------------------------------------------------------------
int JANPACKAMGSolver::getCoarsening(MultigridCoarsening cscheme) const 
{
  switch (cscheme)
  {
  case rs:
    return JP_AMG_RS;
  case cljp:
    return JP_AMG_CLJP;
  case pmis:
    return JP_AMG_PMIS;
  case hmis:
    return JP_AMG_HMIS;
  default:
    warning("Requested Multigrid coarsening unknown. Using Rugen-Stueben.");
    return JP_AMG_RS;
  }
}
//-----------------------------------------------------------------------------
#endif
