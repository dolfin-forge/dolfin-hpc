// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PETSC_KRYLOV_MATRIX_H
#define __DOLFIN_PETSC_KRYLOV_MATRIX_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include <petscmat.h>

#include <dolfin/common/types.h>
#include <dolfin/log/dolfin_log.h>
#include "PETScObject.h"

namespace dolfin
{

  class PETScVector;
  
  /// This class represents a matrix-free matrix of dimension m x m.
  /// It is a simple wrapper for a PETSc shell matrix. The interface
  /// is intentionally simple. For advanced usage, access the PETSc
  /// Mat pointer using the function mat() and use the standard PETSc
  /// interface.
  ///
  /// The class PETScKrylovMatrix enables the use of Krylov subspace
  /// methods for linear systems Ax = b, without having to explicitly
  /// store the matrix A. All that is needed is that the user-defined
  /// PETScKrylovMatrix implements multiplication with vectors. Note that
  /// the multiplication operator needs to be defined in terms of
  /// PETSc data structures (Vec), since it will be called from PETSc.

  class PETScKrylovMatrix : public PETScObject
  {
  public:

    /// Constructor
    PETScKrylovMatrix() = default;

    /// Create a virtual matrix matching the given vectors
    PETScKrylovMatrix(const PETScVector& x, const PETScVector& y);

    /// Destructor
    virtual ~PETScKrylovMatrix();

    /// Initialize virtual matrix matching the given vectors
    void init(const PETScVector& x, const PETScVector& y);

    /// Initialize virtual matrix of given size, all on one processor
    void init(int M, int N);

    /// Return number of rows (dim = 0) or columns (dim = 1) along dimension dim
    auto size(uint dim) const -> uint;

    /// Return PETSc Mat pointer
    auto mat() const -> Mat;

    /// Compute product y = Ax
    virtual void mult(const PETScVector& x, PETScVector& y) const = 0;

    /// Display matrix (sparse output is default)
    void disp(bool sparse = true, int precision = 2) const;

    /// Output
    friend auto operator<< (LogStream& stream, const PETScKrylovMatrix& A) -> LogStream&;

  private:

    // PETSc Mat pointer
    Mat A{ nullptr };

  };

  auto operator<< (LogStream& stream, const PETScKrylovMatrix& A) -> LogStream&;
  
}

#endif

#endif
