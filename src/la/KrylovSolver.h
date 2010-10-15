// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Ola Skavhaug, 2008.
// Modified by Anders Logg, 2008.
//
// First added:  2007-07-03
// Last changed: 2008-06-13

#ifndef __KRYLOV_SOLVER_H
#define __KRYLOV_SOLVER_H

#include <dolfin/config/dolfin_config.h>
#include <dolfin/parameter/Parametrized.h>
#include <dolfin/common/Timer.h>
#include "GenericMatrix.h"
#include "GenericVector.h"
#include "PETScMatrix.h"
#include "PETScVector.h"
#include "JANPACKMat.h"
#include "JANPACKVec.h"
#include "PETScKrylovSolver.h"
#include "EpetraKrylovSolver.h"
#include "JANPACKKrylovSolver.h"
#include "SolverType.h"
#include "PreconditionerType.h"

#ifndef NO_UBLAS
#include "uBlasKrylovSolver.h"
#include "uBlasSparseMatrix.h"
#include "uBlasDenseMatrix.h"
#include "uBlasVector.h"
#endif

namespace dolfin
{

  /// This class defines an interface for a Krylov solver. The underlying 
  /// Krylov solver type is defined in default_type.h.
  
  class KrylovSolver : public Parametrized
  {
  public:
    
    /// Create Krylov solver
    KrylovSolver(SolverType solver_type=default_solver, PreconditionerType pc_type=default_pc)
      : solver_type(solver_type), pc_type(pc_type), ublas_solver(0), petsc_solver(0), 
      epetra_solver(0), janpack_solver(0) {}
    
    /// Destructor
    ~KrylovSolver()
    {
      delete ublas_solver; 
      delete petsc_solver; 
      delete epetra_solver; 
      delete janpack_solver;
    }
    
    /// Solve linear system Ax = b
    uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    { 
      Timer timer("Krylov solver");
#ifndef NO_UBLAS
      if (A.has_type<uBlasSparseMatrix>())
      {
        if (!ublas_solver)
        {
          ublas_solver = new uBlasKrylovSolver(solver_type, pc_type);
          ublas_solver->set("parent", *this);
        }
        return ublas_solver->solve(A.down_cast<uBlasSparseMatrix>(), x.down_cast<uBlasVector>(), b.down_cast<uBlasVector>());
      }

      if (A.has_type<uBlasDenseMatrix>())
      {
        if (!ublas_solver)
        {
          ublas_solver = new uBlasKrylovSolver(solver_type, pc_type);
          ublas_solver->set("parent", *this);
        }
        return ublas_solver->solve(A.down_cast<uBlasDenseMatrix>(), x.down_cast<uBlasVector>(), b.down_cast<uBlasVector>());
      }
#endif
#ifdef HAVE_PETSC
      if (A.has_type<PETScMatrix>())
      {
        if (!petsc_solver)
        {
          petsc_solver = new PETScKrylovSolver(solver_type, pc_type);
          petsc_solver->set("parent", *this);
        }
        return petsc_solver->solve(A.down_cast<PETScMatrix >(), x.down_cast<PETScVector>(), b.down_cast<PETScVector>());
      }
#endif
#ifdef HAVE_TRILINOS
      if (A.has_type<EpetraMatrix>())
      {
        if (!epetra_solver)
        {
          epetra_solver = new EpetraKrylovSolver(solver_type, pc_type);
          epetra_solver->set("parent", *this);
        }
        return epetra_solver->solve(A.down_cast<EpetraMatrix >(), x.down_cast<EpetraVector>(), b.down_cast<EpetraVector>());
      }
#endif
#ifdef HAVE_JANPACK
      if (A.has_type<JANPACKMat>())
      {
	if (!janpack_solver)
	{
	  janpack_solver = new JANPACKKrylovSolver(solver_type, pc_type);
	  janpack_solver->set("parent", *this);
	}
	return janpack_solver->solve(A.down_cast<JANPACKMat>(), x.down_cast<JANPACKVec>(), b.down_cast<JANPACKVec>());
      }
#endif      
      error("No default LU solver for given backend");
      return 0;
    }
    
  private:
    
    // Krylov method
    SolverType solver_type;
    
    // Preconditioner type
    PreconditionerType pc_type;
    
#ifndef NO_UBLAS
    // uBLAS solver
    uBlasKrylovSolver* ublas_solver;
#else
    int* ublas_solver;
#endif

    // PETSc solver
#ifdef HAVE_PETSC
    PETScKrylovSolver* petsc_solver;
#else
    int* petsc_solver;
#endif
#ifdef HAVE_TRILINOS
    EpetraKrylovSolver* epetra_solver;
#else
    int* epetra_solver;
#endif
#ifdef HAVE_JANPACK
    JANPACKKrylovSolver* janpack_solver;
#else
    int* janpack_solver;
#endif
  };
}

#endif
