// Copyright (C) 2006-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MATRIX_H
#define __DOLFIN_MATRIX_H

#include <dolfin/common/Variable.h>
#include <dolfin/main/PE.h>
#include "DefaultFactory.h"
#include "GenericMatrix.h"

#include <sstream>
#include <fstream>

namespace dolfin
{

  /// This class provides the default DOLFIN matrix class,
  /// based on the default DOLFIN linear algebra backend.

  class Matrix : public GenericMatrix, public Variable
  {
  public:

    /// Create empty matrix
    Matrix() : Variable("A", "DOLFIN matrix"), matrix(0)
    { DefaultFactory factory; matrix = factory.createMatrix(); }

    /// Create M x N matrix distributed by default
    Matrix(uint M, uint N) : Variable("A", "DOLFIN matrix"), matrix(0)
    { DefaultFactory factory; matrix = factory.createMatrix(); matrix->init(M, N); }

    /// Create M x N matrix distributed if specified
    Matrix(uint M, uint N, bool distributed) : Variable("A", "DOLFIN matrix"), matrix(0)
    { DefaultFactory factory; matrix = factory.createMatrix(); matrix->init(M, N, distributed); }

    /// Copy constructor
    explicit Matrix(const Matrix& A) : Variable("A", "DOLFIN matrix"),
                                       matrix(A.matrix->copy())
    {}

    /// Destructor
    ~Matrix()
    { delete matrix; }

    //--- Implementation of the GenericTensor interface ---

    /// Initialize zero tensor using sparsity pattern
    /// The tensor is distributed if the sparsity pattern is distributed
    void init(const GenericSparsityPattern& sparsity_pattern)
    { matrix->init(sparsity_pattern); }

    /// Return copy of tensor
    Matrix* copy() const
    { Matrix* A = new Matrix(); delete A->matrix; A->matrix = matrix->copy(); return A; }

    /// Return size of given dimension
    uint size(uint dim) const
    { return matrix->size(dim); }

    /// Set all entries to zero and keep any sparse structure
    void zero()
    { matrix->zero(); }

    /// Finalize assembly of tensor
    void apply(FinalizeType finaltype=FINALIZE)
    { matrix->apply(finaltype); }

    /// Display tensor
    void disp(uint precision=2) const
    { matrix->disp(precision); }

    //--- Implementation of the GenericMatrix interface ---

    /// Initialize M x N matrix and distribute by default
    void init(uint M, uint N)
    { matrix->init(M, N); }

    /// Initialize M x N matrix and distribute if specified
    void init(uint M, uint N, bool distributed)
    { matrix->init(M, N, distributed); }

    /// Get block of values
    void get(real* block, uint m, const uint* rows, uint n, const uint* cols) const
    { matrix->get(block, m, rows, n, cols); }

    /// Set block of values
    void set(const real* block, uint m, const uint* rows, uint n, const uint* cols)
    { matrix->set(block, m, rows, n, cols); }

    /// Add block of values
    void add(const real* block, uint m, const uint* rows, uint n, const uint* cols)
    { matrix->add(block, m, rows, n, cols); }

    /// Return norm of matrix
    double norm(std::string norm_type = "frobenius") const
    { return matrix->norm(norm_type); }

    /// Get non-zero values of given row
    void getrow(uint row, Array<uint>& columns, Array<real>& values) const
    { matrix->getrow(row, columns, values); }

    /// Set values for given row
    void setrow(uint row, const Array<uint>& columns, const Array<real>& values)
    { matrix->setrow(row, columns, values); }

    /// Set given rows to zero
    void zero(uint m, const uint* rows)
    { matrix->zero(m, rows); }

    /// Set given rows to identity matrix
    void ident(uint m, const uint* rows)
    { matrix->ident(m, rows); }

    // Matrix-vector product, y = Ax
    void mult(const GenericVector& x, GenericVector& y, bool transposed=false) const
    { matrix->mult(x, y, transposed); }

    /// Multiply matrix by given number
    const Matrix& operator*= (real a)
    { *matrix *= a; return *this; }

    /// Divide matrix by given number
    const Matrix& operator/= (real a)
    { *matrix /= a; return *this; }

    /// Assignment operator
    const GenericMatrix& operator= (const GenericMatrix& A)
    { *matrix = A; return *this; }

    /// Get number of non-zeros in the matrix
    uint nz() const
    {return  matrix->nz();}

    //--- Special functions ---

    /// Return linear algebra backend factory
    LinearAlgebraFactory& factory() const
    { return matrix->factory(); }

    //--- Special functions, intended for library use only ---

    /// Return concrete instance / unwrap (const)
    const GenericMatrix* instance() const
    { return matrix; }

    /// Return concrete instance / unwrap (non-const version)
    GenericMatrix* instance()
    { return matrix; }

    //--- Special Matrix functions ---

    /// Assignment operator
    const Matrix& operator= (const Matrix& A)
    { *matrix = *A.matrix; return *this; }

    void spy() const;

  private:

    // Pointer to concrete implementation
    GenericMatrix* matrix;

  };

  //--- INLINES --------------------------------------------------------------
  inline void Matrix::spy() const
  {
    if(PE::size() == 1)
    {
      std::stringstream ss;
      ss << "A" << matrix->size(0)*matrix->size(1) << ".xpm" << std::ends;
      std::ofstream Afile(ss.str().c_str());
      real val = 0.;
      std::stringstream xpm;
      xpm << "/* XPM */\n" <<
      "static char * map_xpm = {\n" <<
      "/* width height number_of_colors chars_per_pixel */\n" <<
      "\""<< matrix->size(0) << " "<< matrix->size(1)<<" 2 1\",\n" <<
      "/* intensity levels */\n"<<
      "\"0 c #ffffff\",\n" <<
      "\"1 c none\",\n" <<
      "/* map */\n";
      Afile << xpm.str();
      for (uint i = 0; i < matrix->size(0); ++i)
      {
        std::stringstream row;
        row << "\"";
        for (uint j = 0; j < matrix->size(1); ++j)
        {
          val = matrix->getitem(std::pair<uint,uint>(i,j));
          if(val > 1.0e-14)
          {
            row<< "1";
          }
          else
          {
            row<< "0";
          }
        }
        row << "\",\n";
        Afile << row.str();
      }
      Afile << "};";
    }
  }

}

#endif
