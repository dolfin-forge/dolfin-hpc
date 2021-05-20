// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PETSC_LU_SOLVER_H
#define __DOLFIN_PETSC_LU_SOLVER_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include <dolfin/la/petsc/PETScVector.h>

#include <petscksp.h>
#include <petscmat.h>

namespace dolfin
{

/// Forward declarations
class PETScManager;
class PETScMatrix;
class PETScKrylovMatrix;

//-----------------------------------------------------------------------------

/// This class implements the direct solution (LU factorization) for
/// linear systems of the form Ax = b. It is a wrapper for the LU
/// solver of PETSc.

class PETScLUSolver
{
public:
  /// Constructor
  PETScLUSolver();

  /// Destructor
  ~PETScLUSolver();

  /// Solve linear system Ax = b
  auto solve( const PETScMatrix & A, PETScVector & x, const PETScVector & b )
    -> size_t;

  /// Solve linear system Ax = b
  auto solve( const PETScKrylovMatrix & A,
              PETScVector &             x,
              const PETScVector &       b ) -> size_t;

  /// Display LU solver data
  void disp() const;

private:
  // Create dense copy of virtual matrix
  auto copyToDense( const PETScKrylovMatrix & A ) -> real;

  KSP ksp { nullptr };

  Mat   B { nullptr };
  int * idxm { nullptr };
  int * idxn { nullptr };

  PETScVector e;
  PETScVector y;
};

//-----------------------------------------------------------------------------

} // namespace dolfin

#endif

#endif
