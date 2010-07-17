// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __JANPACK_MAT_H
#define __JANPACK_MAT_H

#ifdef HAS_JANPACK

#include <mat.h>
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
    virtual ~JANPACKMat();

    //--- Implementation of the GenericTensor interface ---

    /// Initialize zero tensor using sparsity pattern
    virtual void init(const GenericSparsityPattern& sparsity_pattern);

    /// Return copy of tensor
    virtual JANPACKMat* copy() const;

    /// Return size of given dimension
    virtual uint size(uint dim) const;

    /// Set all entries to zero and keep any sparse structure
    virtual void zero();

    /// Finalize assembly of tensor
    virtual void apply(FinalizeType finaltype=FINALIZE);

    /// Display tensor
    virtual void disp(uint precision=2) const;

    //--- Implementation of the GenericMatrix interface ---

    /// Initialize M x N matrix
    virtual void init(uint M, uint N);

    /// Get block of values
    virtual void get(real* block, uint m, const uint* rows, uint n, const uint* cols) const;

    /// Set block of values
    virtual void set(const real* block, uint m, const uint* rows, uint n, const uint* cols);

    /// Add block of values
    virtual void add(const real* block, uint m, const uint* rows, uint n, const uint* cols);

    /// Return norm of matrix
    virtual real norm(std::string norm_type = "frobenius") const;

    /// Get non-zero values of given row
    virtual void getrow(uint row, Array<uint>& columns, Array<real>& values) const;

    /// Set values for given row
    virtual void setrow(uint row, const Array<uint>& columns, const Array<real>& values);

    /// Set given rows to zero
    virtual void zero(uint m, const uint* rows);

    /// Set given rows to identity matrix
    virtual void ident(uint m, const uint* rows);

    /// Duplicate matrix
    void dup(GenericMatrix& A); 

    // Matrix-vector product, y = Ax
    virtual void mult(const GenericVector& x, GenericVector& y, bool transposed=false) const;

    /// Multiply matrix by given number
    virtual const JANPACKMat& operator*= (real a);

    /// Divide matrix by given number
    virtual const JANPACKMat& operator/= (real a);

    /// Assignment operator
    virtual const GenericMatrix& operator= (const GenericMatrix& x);

    //--- Special functions ---

    /// Return linear algebra backend factory
    virtual LinearAlgebraFactory& factory() const;

    //--- Special JANPACK functions ---

    /// Return JANPACK Mat_crs pointer;
    Mat_crs*mat() const;
      

    /// Assignment operator
    const JANPACKMat& operator= (const JANPACKMat& x)
    { error("Not implemented."); return *this; }

  private:

    // JANPACK Matrix pointer
    Mat_crs* A;


    Mat_crs _A;
    
    // True if we don't own the matrix A points to
    bool is_view;
    
  };


}

#endif

#endif
