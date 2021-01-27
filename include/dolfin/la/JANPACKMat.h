// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_JANPACK_MAT_H
#define __DOLFIN_JANPACK_MAT_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_JANPACK

#ifdef HAVE_JANPACK_MPI
#define jp_mat_type jp_mat_t
#else
#define jp_mat_type char
#include <janpack/hybrid.h>
#endif

#include <janpack/mat.h>
#include <dolfin/common/Variable.h>
#include "GenericMatrix.h"

namespace dolfin
{
  class GenericSparsityPattern;

  class JANPACKMat: public GenericMatrix, public Variable
  {
  public:

    /// Create empty matrix
    JANPACKMat();

    /// Create M x N matrix
    JANPACKMat(uint M, uint N);

    /// Copy constuctor
    explicit JANPACKMat(const JANPACKMat& A);

    /// Destructor
    ~JANPACKMat();

    //--- Implementation of the GenericTensor interface ---

    /// Initialize zero tensor using sparsity pattern
    void init(const GenericSparsityPattern& sparsity_pattern);

    /// Return copy of tensor
    JANPACKMat* copy() const;

    /// Return size of given dimension
    uint size(uint dim) const;

    /// Set all entries to zero and keep any sparse structure
    void zero();

    /// Finalize assembly of tensor
    void apply(FinalizeType finaltype=FINALIZE);

    /// Display tensor
    void disp(uint precision=2) const;

    //--- Implementation of the GenericMatrix interface ---

    /// Initialize M x N matrix and distribute by default
    void init(uint M, uint N);

    /// Initialize M x N matrix and distribute if specified
    void init(uint M, uint N, bool distributed);

    /// Get block of values
    void get(real* block, uint m, const uint* rows, uint n, const uint* cols) const;

    /// Set block of values
    void set(const real* block, uint m, const uint* rows, uint n, const uint* cols);

    /// Add block of values
    void add(const real* block, uint m, const uint* rows, uint n, const uint* cols);

    /// Return norm of matrix
    real norm(std::string norm_type = "frobenius") const;

    /// Get non-zero values of given row
    void getrow(uint row, std::vector<uint>& columns, std::vector<real>& values) const;

    /// Set values for given row
    void setrow(uint row, const std::vector<uint>& columns, const std::vector<real>& values);

    /// Set given rows to zero
    void zero(uint m, const uint* rows);

    /// Set given rows to identity matrix
    void ident(uint m, const uint* rows);

    // Matrix-vector product, y = Ax
    void mult(const GenericVector& x, GenericVector& y, bool transposed=false) const;

    /// Multiply matrix by given number
    const JANPACKMat& operator*= (real a);

    /// Divide matrix by given number
    const JANPACKMat& operator/= (real a);

    /// Assignment operator
    const GenericMatrix& operator= (const GenericMatrix& x);

    /// Get number of non-zeros in the matrix
    uint nz() const;
    //--- Special functions ---

    /// Return linear algebra backend factory
    LinearAlgebraFactory& factory() const;

    //--- Special JANPACK functions ---

    /// Return JANPACK jp_mat_t pointer;
    jp_mat_type *mat() const;

    /// Assignment operator
    const JANPACKMat& operator= (const JANPACKMat& x)
    { error("Not implemented."); return *this; }

    /// Duplicate matrix
    void dup(const JANPACKMat& A);

  private:

    // JANPACK Matrix pointer
#ifdef HAVE_JANPACK_MPI
    jp_mat_t AA;
    jp_mat_t *A;
#else
    char A[JP_MAT_SIZE_T];
#endif

  };


}

#endif

#endif
