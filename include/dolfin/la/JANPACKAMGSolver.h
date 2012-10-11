// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __JANPACK_AMG_SOLVER_H
#define __JANPACK_AMG_SOLVER_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_JANPACK
#include <dolfin/common/types.h>
#include <dolfin/parameter/Parametrized.h>

#include <janpack/amg_solver.h>

namespace dolfin 
{
  /// Forward declarations
  class JANPACKMat;
  class JANPACKVec;

  class JANPACKAMGSolver : public Parametrized
  {
  public:

    /// Create AMG solver for a particular method and preconditioner
    JANPACKAMGSolver();

    ~JANPACKAMGSolver();
 
    /// Solve linear system Ax = b and return number of iterations
    uint solve(const JANPACKMat& A, JANPACKVec& x, const JANPACKVec& b); 
    
  private:

    /// Get JANPACK Multigrid scheme id
    //    int getType(SolverType method) const;
    

    
  };
}

#endif

#endif
