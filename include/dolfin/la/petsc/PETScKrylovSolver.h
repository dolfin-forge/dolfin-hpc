// Copyright (C) 2004-2005 Johan Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PETSC_KRYLOV_SOLVER_H
#define __DOLFIN_PETSC_KRYLOV_SOLVER_H

#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/la/PreconditionerType.h>
#include <dolfin/la/SolverType.h>

#ifdef HAVE_PETSC

#include <dolfin/la/petsc/PETScPreconditioner.h>

namespace dolfin
{

/// Forward declarations
class PETScMatrix;
class PETScVector;
class PETScKrylovMatrix;

//-----------------------------------------------------------------------------

/// This class implements Krylov methods for linear systems
/// of the form Ax = b. It is a wrapper for the Krylov solvers
/// of PETSc.

class PETScKrylovSolver
{
public:
  /// Create Krylov solver for a particular method and preconditioner
  PETScKrylovSolver( SolverType         method = default_solver,
                     PreconditionerType pc     = PreconditionerType::default_pc );

  /// Create Krylov solver for a particular method and PETScPreconditioner
  PETScKrylovSolver( SolverType            method,
                     PETScPreconditioner & PETScPreconditioner );

  /// Destructor
  ~PETScKrylovSolver();

  /// Solve linear system Ax = b and return number of iterations
  auto solve( const PETScMatrix & A, PETScVector & x, const PETScVector & b )
    -> size_t;

  /// Solve linear system Ax = b and return number of iterations
  auto solve( const PETScKrylovMatrix & A,
              PETScVector &             x,
              const PETScVector &       b ) -> size_t;

  /// Display solver data
  void disp() const;

private:
  /// Initialize KSP solver
  void init( size_t M, size_t N );

  /// Read parameters from database
  void readParameters();

  /// Set solver
  void setSolver();

  /// Set PETScPreconditioner
  void setPETScPreconditioner();

  /// Report the number of iterations
  void writeReport( int num_iterations );

  /// Get PETSc method identifier
  auto getType( SolverType method ) const -> KSPType;

  /// Krylov method
  SolverType method;

  /// PETSc PETScPreconditioner
  PreconditionerType pc_petsc;

  /// DOLFIN PETScPreconditioner
  PETScPreconditioner * pc_dolfin { nullptr };

  /// PETSc solver pointer
  KSP ksp { nullptr };

  /// Size of old system (need to reinitialize when changing)
  size_t M_ { 0 };
  size_t N_ { 0 };

  /// True if we have read parameters
  bool parameters_read { false };

  // FIXME: Required to avoid PETSc bug with Hypre. See explanation inside
  //        PETScKrylovSolver:init(). Can be removed when PETSc is patched.
  bool pc_set { false };
};

//-----------------------------------------------------------------------------

} // namespace dolfin

#endif

#endif
