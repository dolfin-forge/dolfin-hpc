// Copyright (C) 2005 Johan Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PETSC_PRECONDITIONER_H
#define __DOLFIN_PETSC_PRECONDITIONER_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include <dolfin/la/petsc/PETScObject.h>

#include <petscksp.h>
#include <petscpc.h>

namespace dolfin
{

class PETScVector;
enum class PreconditionerType;

/// This class specifies the interface for user-defined Krylov
/// method PETScPreconditioners. A user wishing to implement her own
/// PETScPreconditioner needs only supply a function that approximately
/// solves the linear system given a right-hand side.

class PETScPreconditioner : public PETScObject
{
public:
  /// Constructor
  PETScPreconditioner() = default;

  /// Destructor
  virtual ~PETScPreconditioner() = default;

  static void setup( const KSP ksp, PETScPreconditioner & pc );

  /// Solve linear system approximately for given right-hand side b
  virtual void solve( PETScVector & x, const PETScVector & b ) = 0;

  /// Friends
  friend class PETScKrylovSolver;

protected:
  PC petscpc;

private:
  static auto PCApply( PC pc, Vec x, Vec y ) -> int;
  static auto PCCreate( PC pc ) -> int;

  /// Return PETSc PETScPreconditioner type
  static auto getType( PreconditionerType pc ) -> PCType;
};

}

#endif

#endif
