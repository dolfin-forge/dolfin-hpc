// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TRILINOS_KRYLOV_SOLVER_H
#define __DOLFIN_TRILINOS_KRYLOV_SOLVER_H

#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/la/PreconditionerType.h>
#include <dolfin/la/SolverType.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosPreconditioner.h>
#include <dolfin/la/trilinos/TrilinosVector.h>

#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <Ifpack2_Factory.hpp>

namespace dolfin
{

namespace trilinos
{

/// Forward declarations
class Matrix;
class Vector;

//-----------------------------------------------------------------------------

/// This class implements Krylov methods for linear systems of the form Ax = b.
/// It is a wrapper for the Krylov solvers of trilinos

class KrylovSolver
{
public:
  /// local index type
  using LO = int;

  /// global index type
  using GO = size_t;

  /// Tpetra operator type
  using TPOperator = Tpetra::Operator< real, LO, GO, Vector::TPVector::node_type >;

  /// Belos problem type
  using BLSProblem = Belos::LinearProblem< real, Vector::TPVector, TPOperator >;

  /// Belos Solver RCP
  using BLSSolver  = Belos::SolverManager<real, Vector::TPVector, TPOperator >;

public:
  /// Create Krylov solver for a particular method and preconditioner
  KrylovSolver( SolverType         method = default_solver,
                PreconditionerType pc     = PreconditionerType::default_pc );

  /// Create Krylov solver for a particular method and trilinos::Preconditioner
  KrylovSolver( SolverType                 type,
                trilinos::Preconditioner & Preconditioner );

  /// Destructor
  ~KrylovSolver();

  /// Solve linear system Ax = b and return number of iterations
  auto solve( trilinos::Matrix const & A,
              trilinos::Vector &       x,
              trilinos::Vector const & b ) -> size_t;

  /// Display solver data
  void disp() const;

private:
  friend trilinos::Preconditioner;

  /// Initialize KSP solver
  void init();

  /// Read parameters from database
  void readParameters();

  /// Set solver
  void setSolver();

  /// Set PETScPreconditioner
  void setPETScPreconditioner();

  /// Report the number of iterations
  void writeReport( int num_iterations );

  /// Get PETSc method identifier
  // auto getType( SolverType method ) const -> KSPType;

  /// Krylov method
  SolverType solver_type_;

  // Belos solver pointer
  Teuchos::RCP< BLSSolver > solver_;

  /// Preconditioner Type
  PreconditionerType pc_type_;

  // The preconditioner, if any
  Preconditioner * pc_ { nullptr };

  // Container for the problem
  Teuchos::RCP< BLSProblem > problem_;
};

//-----------------------------------------------------------------------------

} // namespace trilinos

} // namespace dolfin

#endif // HAVE_TRILINOS

#endif // __DOLFIN_TRILINOS_KRYLOV_SOLVER_H
