// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __AMG_SOLVER_H
#define __AMG_SOLVER_H

#include <dolfin/config/dolfin_config.h>
#include <dolfin/parameter/Parametrized.h>
#include <dolfin/common/Timer.h>
#include "GenericMatrix.h"
#include "GenericVector.h"
#include "JANPACKMat.h"
#include "JANPACKVec.h"
#include "JANPACKAMGSolver.h"
#include "MultigridScheme.h"

namespace dolfin
{

  /// This class defines an interface for a AMG solver. 
  
  class AMGSolver : public Parametrized
  {
  public:
    
    /// Create Krylov solver
    AMGSolver(MultigridScheme scheme_type=default_scheme) 
      : scheme_type(scheme_type), janpack_solver(0) {}
    
    /// Destructor
    ~AMGSolver()
    {
      delete janpack_solver;
    }
    
    /// Solve linear system Ax = b
    uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    { 
      Timer timer("AMG solver");

#ifdef HAVE_JANPACK
      if (A.has_type<JANPACKMat>())
      {
	if (!janpack_solver)
	{
	  janpack_solver = new JANPACKAMGSolver(scheme_type);
	  janpack_solver->set("parent", *this);
	}
	return janpack_solver->solve(A.down_cast<JANPACKMat>(), x.down_cast<JANPACKVec>(), b.down_cast<JANPACKVec>());
      }
#endif      
      error("No AMG solver for given backend");
      return 0;
    }
    
  private:

    // Multigrid Scheme
    MultigridScheme scheme_type;


#ifdef HAVE_JANPACK
    JANPACKAMGSolver* janpack_solver;
#else
    int* janpack_solver;
#endif
  };
}

#endif
