// Copyright (C) 2004-2008 Johan Hoffman, Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PETSC_MATRIX_H
#define __DOLFIN_PETSC_MATRIX_H

#include <dolfin/config/dolfin_config.h>

#include <dolfin/common/Array.h>
#include <dolfin/common/Variable.h>
#include <dolfin/la/GenericMatrix.h>

#ifdef HAVE_PETSC

#include <dolfin/la/PETScObject.h>

#include <petscmat.h>

namespace dolfin
{

class PETScVector;

/// This class provides a simple matrix class based on PETSc.
/// It is a simple wrapper for a PETSc matrix pointer (Mat)
/// implementing the GenericMatrix interface.
///
/// The interface is intentionally simple. For advanced usage,
/// access the PETSc Mat pointer using the function mat() and
/// use the standard PETSc interface.

class PETScMatrix : public GenericMatrix, public PETScObject, public Variable
{

public:


  enum Norm
  {
    l1, linf, frobenius
  };

  /// Create empty matrix
  explicit PETScMatrix();

  /// Create matrix of local dimension M x N
  PETScMatrix(uint M, uint N, bool distributed = true);

  /// Copy constructor
  explicit PETScMatrix(const PETScMatrix& A);

  /// Create matrix from given PETSc Mat pointer
  explicit PETScMatrix(Mat A);

  /// Destructor
  ~PETScMatrix();

  //--- Implementation of the GenericTensor interface ---

  /// Initialize zero tensor using sparsity pattern
  void init(const GenericSparsityPattern& sparsity_pattern);

  /// Return copy of tensor
  PETScMatrix* copy() const;

  /// Return size of given dimension
  uint size(uint dim) const;

  /// Set all entries to zero and keep any sparse structure
  void zero();

  /// Finalize assembly of tensor
  void apply(FinalizeType final = FINALIZE);

  /// Display tensor
  void disp(uint precision = 0) const;

  //--- Implementation of the GenericMatrix interface --

  /// Initialize matrix of local dimension M x N, distributed by default
  void init(uint M, uint N);

  /// Initialize matrix of local dimension M x N, distributed if specified
  void init(uint M, uint N, bool distributed);

  /// Get block of values
  void get(real* block, uint m, const uint* rows, uint n,
           const uint* cols) const;

  /// Set block of values
  void set(const real* block, uint m, const uint* rows, uint n,
           const uint* cols);

  /// Add block of values
  void add(const real* block, uint m, const uint* rows, uint n,
           const uint* cols);

  /// Return norm of matrix
  real norm(std::string norm_type = "frobenius") const;

  /// Get non-zero values of given row
  void getrow(uint row, Array<uint>& columns, Array<real>& values) const;

  /// Set values for given row
  void setrow(uint row, const Array<uint>& columns, const Array<real>& values);

  /// Set given rows to zero
  void zero(uint m, const uint* rows);

  /// Set given rows to identity matrix
  void ident(uint m, const uint* rows);

  /// Duplicate matrix
  void dup(GenericMatrix& A);

  // Matrix-vector product, y = Ax
  void mult(const GenericVector& x, GenericVector& y,
            bool transposed = false) const;

  /// Multiply matrix by given number
  const PETScMatrix& operator*=(real a);

  /// Divide matrix by given number
  const PETScMatrix& operator/=(real a);

  /// Assignment operator
  const GenericMatrix& operator=(const GenericMatrix& A);

  /// Get number of non-zeros in the matrix
  uint nz() const;

  //--- Special functions ---

  /// Return linear algebra backend factory
  LinearAlgebraFactory& factory() const;

  //--- Special PETScFunctions ---

  /// Return PETSc Mat pointer
  Mat mat() const;

  /// Return norm of matrix
  real norm(const Norm type = l1) const;

  /// Assignment operator
  const PETScMatrix& operator=(const PETScMatrix& A);

  /// Matrix axpy, Y = a X+ Y
  const PETScMatrix& operator+=(const PETScMatrix& A);

private:

  //
  void clear();

  // Initialize M x N matrix with a given number of nonzeros per row
  void init(uint M, uint N, const uint* nz);

  // Initialize M x N matrix with a given number of nonzeros per row diagonal
  // and off-diagonal
  void init(uint M, uint N, const uint* d_nzrow, const uint* o_nzrow);

  ///
  void getrows_offproc(_ordered_set<uint> const& rows);

  // Print info
  void print(MatInfo const& info) const;

  // Matrix
  Mat A;

  // Sub-matrices
  Mat * AA_sub;

  // True if the matrix is distributed
  bool is_distributed_;

  PetscInt rstart_;
  PetscInt rend_;

  _map<int, int> mapping_;

};

//-----------------------------------------------------------------------------
inline void PETScMatrix::init( uint M, uint N )
{
  init( M, N, true );
}
//-----------------------------------------------------------------------------
inline uint PETScMatrix::size( uint dim ) const
{
  int M = 0;
  int N = 0;
  MatGetSize( A, &M, &N );
  return ( dim == 0 ? M : N );
}
//-----------------------------------------------------------------------------
inline uint PETScMatrix::nz() const
{
  MatInfo info;
  MatGetInfo( A, MAT_GLOBAL_SUM, &info );

  return info.nz_used;
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::get( real* block,
                              uint m, uint const* rows,
                              uint n, uint const* cols ) const
{
  dolfin_assert(A);
  MatGetValues(A, static_cast<int>(m),
               reinterpret_cast<int*>(const_cast<uint*>(rows)),
               static_cast<int>(n),
               reinterpret_cast<int*>(const_cast<uint*>(cols)), block);
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::set( real const* block,
                              uint m, uint const* rows,
                              uint n, uint const* cols )
{
  dolfin_assert(A);
  MatSetValues(A, static_cast<int>(m),
               reinterpret_cast<int*>(const_cast<uint*>(rows)),
               static_cast<int>(n),
               reinterpret_cast<int*>(const_cast<uint*>(cols)), block,
               INSERT_VALUES);
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::add( real const* block,
                              uint m, uint const* rows,
                              uint n, uint const* cols )
{
  dolfin_assert(A);
  MatSetValues(A, static_cast<int>(m),
               reinterpret_cast<int*>(const_cast<uint*>(rows)),
               static_cast<int>(n),
               reinterpret_cast<int*>(const_cast<uint*>(cols)), block,
               ADD_VALUES);
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::setrow( uint row, const Array<uint> & columns,
                                 const Array<real>& values )
{
  set( values.data(), 1, &row, columns.size(), columns.data() );
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::apply( FinalizeType finaltype )
{
  if ( finaltype == FINALIZE )
  {
    MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );
  }
  else if ( finaltype == FLUSH )
  {
    MatAssemblyBegin( A, MAT_FLUSH_ASSEMBLY );
    MatAssemblyEnd( A, MAT_FLUSH_ASSEMBLY );
  }
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::zero()
{
  MatZeroEntries( A );
}
//-----------------------------------------------------------------------------
inline const PETScMatrix & PETScMatrix::operator+=( const PETScMatrix & A )
{
  dolfin_assert( this->A );
  MatAXPY( this->A, 1.0, A.A, SAME_NONZERO_PATTERN );
  return *this;
}
//-----------------------------------------------------------------------------
inline const PETScMatrix & PETScMatrix::operator*=( real a )
{
  dolfin_assert( A );
  MatScale( A, a );
  return *this;
}
//-----------------------------------------------------------------------------
inline const PETScMatrix & PETScMatrix::operator/=( real a )
{
  dolfin_assert( A );
  MatScale( A, 1.0 / a );
  return *this;
}
//-----------------------------------------------------------------------------
inline const GenericMatrix & PETScMatrix::operator=( const GenericMatrix & A )
{
  if ( &A != this )
  {
    MatCopy( A.down_cast< PETScMatrix >().A, this->A, SAME_NONZERO_PATTERN );
  }
  return *this;
}
//-----------------------------------------------------------------------------
inline const PETScMatrix & PETScMatrix::operator=( const PETScMatrix & A )
{
  if ( &A != this )
  {
    MatCopy( A.A, ( this->A ), SAME_NONZERO_PATTERN );
  }
  return *this;
}
//-----------------------------------------------------------------------------
inline void PETScMatrix::dup( GenericMatrix & A )
{
  MatDuplicate( A.down_cast< PETScMatrix >().A, MAT_COPY_VALUES, &this->A );
}
//-----------------------------------------------------------------------------
inline Mat PETScMatrix::mat() const
{
  return A;
}

} /* namespace dolfin */

#endif /* HAVE_PETSC */

#endif /* __DOLFIN_PETSC_MATRIX_H */
