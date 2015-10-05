// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Ola Skavhaug, 2008.
// Modified by Anders Logg, 2008.
// Modified by Niclas Jansson, 2010-2011.
//
// First added:  2007-07-03
// Last changed: 2011-06-10

#ifndef __DOLFIN_KRYLOV_SOLVER_H
#define __DOLFIN_KRYLOV_SOLVER_H

#include <dolfin/config/dolfin_config.h>

#include <dolfin/la/SolverType.h>
#include <dolfin/la/PreconditionerType.h>
#include <dolfin/parameter/Parametrized.h>

#include "PETScMatrix.h"
#include "PETScVector.h"
#include "JANPACKMat.h"
#include "JANPACKVec.h"
#include "PETScKrylovSolver.h"
#include "JANPACKKrylovSolver.h"

namespace dolfin
{

class GenericMatrix;
class GenericVector;

class JANPACKKrylovSolver;
class JANPACKMat;
class JANPACKVec;
class PETScKrylovSolver;
class PETScMatrix;
class PETScVector;

/// This class defines an interface for a Krylov solver. The underlying
/// Krylov solver type is defined in default_type.h.

class KrylovSolver : public Parametrized
{

public:

  /// Create Krylov solver
  KrylovSolver(SolverType solver_type = default_solver,
               PreconditionerType pc_type = default_pc);

  /// Destructor
  ~KrylovSolver();

  /// Solve linear system Ax = b
  uint solve(GenericMatrix const& A, GenericVector& x, GenericVector const& b);

private:

  /// Krylov method
  SolverType solver_type;

  /// Preconditioner type
  PreconditionerType pc_type;

  /// PETSc solver
#ifdef HAVE_PETSC
  PETScKrylovSolver * petsc_solver;
#else
  int * petsc_solver;
#endif
#ifdef HAVE_JANPACK
  JANPACKKrylovSolver * janpack_solver;
#else
  int * janpack_solver;
#endif
};

}

#endif
