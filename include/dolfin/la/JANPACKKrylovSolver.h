// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __DOLFIN_JANPACK_KRYLOV_SOLVER_H
#define __DOLFIN_JANPACK_KRYLOV_SOLVER_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_JANPACK
#include <dolfin/common/types.h>
#include <dolfin/parameter/Parametrized.h>
#include "SolverType.h"
#include "PreconditionerType.h"


#ifndef HAVE_JANPACK_MPI
#include <janpack/hybrid.h>
#endif

#include <janpack/krylov_solver.h>

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

    ~JANPACKKrylovSolver();
 
    /// Solve linear system Ax = b and return number of iterations
    uint solve(const JANPACKMat& A, JANPACKVec& x, const JANPACKVec& b); 
    
  private:

    /// Krylov method
    SolverType method;

    /// Get JANPACK krylov method id
    int getType(SolverType method) const;

    /// DOLFIN PETScPreconditioner
    PreconditionerType pc_janpack;
    
    bool ksp_init;

    // JANPACK ksp data
#ifdef HAVE_JANPACK_MPI
    jp_ksp_t _ksp;
    jp_ksp_t *ksp;
#else
    char ksp[JP_KSP_SIZE_T];
#endif
    
  };
}

#endif

#endif
