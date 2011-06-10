// Copyright (C) 2007-2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Ola Skavhaug, 2008.
// Modified by Dag Lindbo, 2008.
// Modified by Anders Logg, 2008.
// Modified by Kent-Andre Mardal, 2008.
// Modified by Niclas Jansson, 2010-2011.
//
// First added:  2007-07-03
// Last changed: 2011-06-10

#ifndef __LU_SOLVER_H
#define __LU_SOLVER_H

#include <dolfin/config/dolfin_config.h>
#include <dolfin/parameter/Parametrized.h>
#include <dolfin/common/Timer.h>
#include "GenericMatrix.h"
#include "GenericVector.h"
#include "PETScLUSolver.h"
#include "PETScMatrix.h"


namespace dolfin
{

  class LUSolver : public Parametrized
  {

  /// LU solver for the built-in LA backends. 
    
  public:

    LUSolver() : petsc_solver(0) {}
    
    ~LUSolver() 
    { 
      delete petsc_solver; 
    }
    
    uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    {
      Timer timer("LU solver");

#ifdef HAVE_PETSC
      if (A.has_type<PETScMatrix>()) 
      {
        if (!petsc_solver)
        {
          petsc_solver = new PETScLUSolver();
          petsc_solver->set("parent", *this);
        }
        return petsc_solver->solve(A.down_cast<PETScMatrix>(), x.down_cast<PETScVector>(), b.down_cast<PETScVector>());
      }
#endif
      error("No default LU solver for given backend");
      return 0;
    }

    uint factorize(const GenericMatrix& A)
    {

      error("No matrix factorization for given backend.");
      return 0;

    }
    
    uint factorized_solve(GenericVector& x, const GenericVector& b)
    {
      error("No factorized LU solver for given backend.");
      return 0;
    }

  private:

    // PETSc Solver
#ifdef HAVE_PETSC
    PETScLUSolver* petsc_solver;
#else
    int* petsc_solver;
#endif

  };
}

#endif
