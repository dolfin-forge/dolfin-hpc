// Copyright (C) 2004-2008 Johan Hoffman, Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2005-2007.
// Modified by Andy R. Terrel 2005.
// Modified by Ola Skavhaug 2007.
// Modified by Magnus Vikstrøm 2007-2008.
// Modified by Niclas Jansson 2008-2010.
//
// First added:  2004
// Last changed: 2010-01-03

#ifdef HAS_PETSC

#include <iostream>
#include <sstream>
#include <iomanip>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/Array.h>
#include "PETScVector.h"
#include "PETScMatrix.h"
#include "GenericSparsityPattern.h"
#include "SparsityPattern.h"
#include "PETScFactory.h"
#include <dolfin/main/MPI.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix(const Type type):
    Variable("A", "a sparse matrix"),
    A(0), is_view(false), _type(type), sub(false), block_size(0)
{
  // Check type
  checkType();
}
//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix(Mat A):
    Variable("A", "a sparse matrix"),
    A(A), is_view(true), _type(default_matrix), block_size(0)
{
  // FIXME: get PETSc matrix type and set
  _type = default_matrix;
}
//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix(uint M, uint N, Type type):
    Variable("A", "a sparse matrix"),
    A(0), is_view(false), _type(type), block_size(0)
{
  // Check type
  checkType();

  // Create PETSc matrix
  init(M, N);
}
//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix(const PETScMatrix& A):
  Variable("A", "PETSc matrix"),
  A(0), is_view(false), _type(A._type), block_size(0)
{
  *this = A;
}
//-----------------------------------------------------------------------------
PETScMatrix::~PETScMatrix()
{
  if (A && !is_view)
    MatDestroy(A);

  if(sub)
    MatDestroy(AA_sub[0]); 

  // FIXME destroy sub
}
//-----------------------------------------------------------------------------
void PETScMatrix::init(uint M, uint N)
{
  // Free previously allocated memory if necessary
  if ( A )
    MatDestroy(A);
  
  // FIXME: maybe 50 should be a parameter?
  // FIXME: it should definitely be a parameter

  // Create a sparse matrix in compressed row format
  if (dolfin::MPI::numProcesses() > 1)
  {
    // Create PETSc parallel matrix with a guess for number of diagonal (50 in this case) 
    // and number of off-diagonal non-zeroes (50 in this case).
    // Note that guessing too high leads to excessive memory usage.
    // In order to not waste any memory one would need to specify d_nnz and o_nnz.
#ifdef HAS_MPI
    MatCreateMPIAIJ(dolfin::MPI::DOLFIN_COMM, PETSC_DECIDE, PETSC_DECIDE, 
		    M, N, 120, PETSC_NULL, 120, PETSC_NULL, &A);
#endif
  }
  else
  {
    // Create PETSc sequential matrix with a guess for number of non-zeroes (50 in thise case)
    MatCreateSeqAIJ(PETSC_COMM_SELF, M, N, 50, PETSC_NULL, &A);

    setType();
    MatSetOption(A, MAT_KEEP_ZEROED_ROWS);
    MatSetFromOptions(A);
  }
}
//-----------------------------------------------------------------------------
void PETScMatrix::init(uint M, uint N, const uint* nz)
{
  // Free previously allocated memory if necessary
  if ( A )
    MatDestroy(A);

  // Create a sparse matrix in compressed row format
  if (dolfin::MPI::numProcesses() > 1)
  {
    // Create PETSc parallel matrix with a guess for number of diagonal (50 in this case) 
    // and number of off-diagonal non-zeroes (50 in this case).
    // Note that guessing too high leads to excessive memory usage.
    // In order to not waste any memory one would need to specify d_nnz and o_nnz.
#ifdef HAS_MPI
    MatCreateMPIAIJ(dolfin::MPI::DOLFIN_COMM, PETSC_DECIDE, PETSC_DECIDE, 
		    M, N, 120, PETSC_NULL, 120, PETSC_NULL, &A);
#endif
    //MatSetFromOptions(A);
    //MatSetOption(A, MAT_KEEP_ZEROED_ROWS);
    //MatZeroEntries(A);
  }
  else
  {
    //    MatCreateSeqBAIJ(MPI_Comm comm,PetscInt bs,PetscInt m,PetscInt n,PetscInt nz,const PetscInt nnz[],Mat *A);
    

    if ( block_size ) 
    {      
      MatCreateSeqBAIJ(PETSC_COMM_SELF, block_size, (int) M, (int) N, 1, PETSC_NULL, &A);
      MatSetOption(A, MAT_KEEP_ZEROED_ROWS);
      MatSetOption(A, MAT_USE_HASH_TABLE);
      MatSetFromOptions(A);
      MatZeroEntries(A);
      //    MatCreateSeqBAIJ(PETSC_COMM_SELF, 1, (int) M, (int) N, 5, PETSC_NULL, &A);
      // Create PETSc sequential matrix with a guess for number of non-zeroes (50 in thise case)
    }
    else 
    {
      MatCreate(PETSC_COMM_SELF, &A);
      MatSetSizes(A,  PETSC_DECIDE,  PETSC_DECIDE, M, N);
      setType();
      MatSetOption(A, MAT_KEEP_ZEROED_ROWS);
      MatSetFromOptions(A);
      MatSeqAIJSetPreallocation(A, PETSC_DEFAULT, (int*)nz);
      MatZeroEntries(A);
      } 
  }
}
//-----------------------------------------------------------------------------
void PETScMatrix::init(uint M, uint N, const uint* d_nzrow, const uint* o_nzrow)
{
  // Free previously allocated memory if necessary
  if ( A )
    MatDestroy(A);
  // Create PETSc parallel matrix with a guess for number of diagonal (50 in this case) 
  // and number of off-diagonal non-zeroes (50 in this case).
  // Note that guessing too high leads to excessive memory usage.
  // In order to not waste any memory one would need to specify d_nnz and o_nnz.

  //  MatCreateMPIAIJ(dolfin::MPI::DOLFIN_COMM, PETSC_DECIDE, PETSC_DECIDE, 
		  //		  M, N, PETSC_NULL, (int*)d_nzrow, PETSC_NULL, (int*)o_nzrow, &A);
#ifdef HAS_MPI
  MatCreateMPIAIJ(MPI::DOLFIN_COMM, M, N, PETSC_DETERMINE, PETSC_DETERMINE, PETSC_NULL, (int*)d_nzrow, PETSC_NULL, (int*)o_nzrow, &A);
#endif
  //MatCreateMPIAIJ(MPI::DOLFIN_COMM, (int) M, (int) N, PETSC_DETERMINE, PETSC_DETERMINE, 90, PETSC_NULL, 90, PETSC_NULL, &A);


  //MatSetOption(A, MAT_COLUMNS_SORTED); // assert("it's going to be ok");
  MatSetOption(A, MAT_KEEP_ZEROED_ROWS);
  MatSetFromOptions(A);
  MatZeroEntries(A);

}
//-----------------------------------------------------------------------------
void PETScMatrix::init(const GenericSparsityPattern& sparsity_pattern)
{

  if (dolfin::MPI::numProcesses() > 1)
  {
    uint p = dolfin::MPI::processNumber();
    const SparsityPattern& spattern = reinterpret_cast<const SparsityPattern&>(sparsity_pattern);
    uint local_size = spattern.numLocalRows(p);
    uint* d_nzrow = new uint[local_size];
    uint* o_nzrow = new uint[local_size];
    spattern.numNonZeroPerRow(p, d_nzrow, o_nzrow);
    uint spattern_size_0 = spattern.size(0);
    uint spattern_size_1 = spattern.size(1);
    const_cast<SparsityPattern&>(spattern).clear();
    init(spattern_size_0, spattern_size_1, d_nzrow, o_nzrow);
    delete [] d_nzrow;
    delete [] o_nzrow;
  }
  else
  {

    const SparsityPattern& spattern = reinterpret_cast<const SparsityPattern&>(sparsity_pattern);
    uint* nzrow = new uint[spattern.size(0)];  
    spattern.numNonZeroPerRow(nzrow);
    init(spattern.size(0), spattern.size(1), nzrow);
    delete [] nzrow;
  }
}
//-----------------------------------------------------------------------------
PETScMatrix* PETScMatrix::copy() const
{

  PETScMatrix* mcopy = new PETScMatrix();
  
  MatDuplicate(A, MAT_COPY_VALUES, &(mcopy->A));

  return mcopy;
}
//-----------------------------------------------------------------------------
dolfin::uint PETScMatrix::size(uint dim) const
{
  int M = 0;
  int N = 0;
  MatGetSize(A, &M, &N);
  return (dim == 0 ? M : N);
}
//-----------------------------------------------------------------------------
void PETScMatrix::get(real* block,
                      uint m, const uint* rows,
                      uint n, const uint* cols) const
{
  dolfin_assert(A);
  MatGetValues(A,
               static_cast<int>(m), reinterpret_cast<int*>(const_cast<uint*>(rows)),
               static_cast<int>(n), reinterpret_cast<int*>(const_cast<uint*>(cols)),
               block);
}
//-----------------------------------------------------------------------------
void PETScMatrix::set(const real* block,
                      uint m, const uint* rows,
                      uint n, const uint* cols)
{
  dolfin_assert(A);



  if (block_size) 
    
    MatSetValuesBlocked(A,
			static_cast<int>(m) / block_size,
			reinterpret_cast<int*>(const_cast<uint*>(rows)),
			static_cast<int>(n) / block_size,
			reinterpret_cast<int*>(const_cast<uint*>(cols)),
			block, INSERT_VALUES);

  else 

    MatSetValues(A,
		 static_cast<int>(m), reinterpret_cast<int*>(const_cast<uint*>(rows)),
		 static_cast<int>(n), reinterpret_cast<int*>(const_cast<uint*>(cols)),
		 block, INSERT_VALUES);

}
//-----------------------------------------------------------------------------
void PETScMatrix::add(const real* block,
                      uint m, const uint* rows,
                      uint n, const uint* cols)
{
  dolfin_assert(A);

  if ( block_size )  
    MatSetValuesBlocked(A,
			static_cast<int>(m) / block_size,
			reinterpret_cast<int*>(const_cast<uint*>(rows)),
			static_cast<int>(n) / block_size,
			reinterpret_cast<int*>(const_cast<uint*>(cols)),
			block, ADD_VALUES);
  
  else
  
    MatSetValues(A,
		 static_cast<int>(m), reinterpret_cast<int*>(const_cast<uint*>(rows)),
		 static_cast<int>(n), reinterpret_cast<int*>(const_cast<uint*>(cols)),
		 block, ADD_VALUES);
}
//-----------------------------------------------------------------------------
real PETScMatrix::norm(std::string norm_type) const
{
  real norm;
  if (norm_type == "l1")
    MatNorm(A, NORM_1, &norm);
  else if (norm_type == "linf")
    MatNorm(A, NORM_INFINITY, &norm);
  else if (norm_type == "frobenius")
    MatNorm(A, NORM_FROBENIUS, &norm);
  else
  {
    error("Unknown norm type in uPETScMatrix.");
    return 0.0;
  }
  return norm;
}
//-----------------------------------------------------------------------------
void PETScMatrix::getrow(uint row,
                         Array<uint>& columns,
                         Array<real>& values) const
{
  const int *cols = 0;
  const double *vals = 0;
  int ncols = 0;

  int m, n;
  MatGetOwnershipRange(A, &m, &n);
    
  // Two different cases, local and non local rows.
  if( row >= (uint) m && row < (uint) n) {

    MatGetRow(A, row, &ncols, &cols, &vals);
    
    // Assign values to Arrays
    columns.assign(reinterpret_cast<uint *>(const_cast<int*>(cols)),
		   reinterpret_cast<uint*>(const_cast<int*>(cols+ncols)));
    values.assign(vals, vals+ncols);
    
    MatRestoreRow(A, row, &ncols, &cols, &vals);

  }
  else {
    if(!sub)
      error("No ghosted Processor rows");
#if sun
    std::map<int, int>::const_iterator it = mapping.find(row);    
#else
    std::map<const int, int>::const_iterator it = mapping.find(row);    
#endif
    MatGetRow(AA_sub[0], it->second, &ncols, &cols, &vals);
    columns.assign(reinterpret_cast<uint*>(const_cast<int*>(cols)), 
		   reinterpret_cast<uint*>(const_cast<int*>(cols+ncols)));
    values.assign(vals, vals+ncols);
    MatRestoreRow(AA_sub[0], it->second, &ncols, &cols, &vals);
  }
  
} 
//-----------------------------------------------------------------------------
void PETScMatrix::getrows_offproc(std::set<uint> rows)
{

  if(dolfin::MPI::numProcesses() == 1)
    return;

  int *_cols = new int[size(0)];
  int *_rows = new int[rows.size()];

  int m, n;
  MatGetOwnershipRange(A, &m, &n);

  std::set<uint>::iterator it;

  mapping.clear();

  uint i = 0;
  for(it = rows.begin(); it != rows.end(); it++) {
    if( !(*it >= (uint) m && *it < (uint) n)) {
      _rows[i] = *it;
      mapping[*it] = i++;
    }
  }
 
  for(uint j = 0; j < size(0); j++)
    _cols[j] = j;

  
  IS irow, icol;
  ISCreateGeneral(PETSC_COMM_SELF, i, _rows, &irow);
  ISCreateGeneral(PETSC_COMM_SELF, size(0), _cols, &icol);

  if(!sub)
    MatGetSubMatrices(A, 1, &irow, &icol, MAT_INITIAL_MATRIX, &AA_sub);    
  else
    MatGetSubMatrices(A, 1, &irow, &icol, MAT_REUSE_MATRIX, &AA_sub);    

  sub = true;
  
  ISDestroy(irow);
  ISDestroy(icol);

  delete[] _cols;
  delete[] _rows;
  
}
//-----------------------------------------------------------------------------
void PETScMatrix::setrow(uint row,
                         const Array<uint>& columns,
                         const Array<real>& values)
{
  // Check size of arrays
  if (columns.size() != values.size())
    error("Number of columns and values don't match for setrow() operation.");

  // Handle case n = 0
  const uint n = columns.size();
  if (n == 0)
    return;

  // Assign values to arrays
  uint* cols = new uint[n];
  real* vals = new real[n];
  for (uint j = 0; j < n; j++)
  {
    cols[j] = columns[j];
    vals[j] = values[j];
  }
  
  // Set values
  set(vals, 1, &row, n, cols);
  
  // Free temporary storage
  delete [] cols;
  delete [] vals;
}
//-----------------------------------------------------------------------------
void PETScMatrix::zero(uint m, const uint* rows)
{
  IS is = 0;
  ISCreateGeneral(PETSC_COMM_SELF, static_cast<int>(m), reinterpret_cast<int*>(const_cast<uint*>(rows)), &is);
  PetscScalar null = 0.0;
  MatZeroRowsIS(A, is, null);
  ISDestroy(is);
}
//-----------------------------------------------------------------------------
void PETScMatrix::ident(uint m, const uint* rows)
{
  IS is = 0;
  ISCreateGeneral(PETSC_COMM_SELF, static_cast<int>(m), reinterpret_cast<int*>(const_cast<uint*>(rows)), &is);
  PetscScalar one = 1.0;
  MatZeroRowsIS(A, is, one);
  ISDestroy(is);
}
//-----------------------------------------------------------------------------
void PETScMatrix::mult(const GenericVector& x, GenericVector& y, bool transposed) const
{
  const PETScVector& xx = x.down_cast<PETScVector>();  
  PETScVector& yy = y.down_cast<PETScVector>();

  if (transposed)
  { 
    if (size(0) != xx.size()) error("Matrix and vector dimensions don't match for matrix-vector product.");
    yy.init(size(1));
    MatMultTranspose(A, xx.vec(), yy.vec());
  }
  else {
    if (size(1) != xx.size()) error("Matrix and vector dimensions don't match for matrix-vector product.");
    yy.init(size(0));
    MatMult(A, xx.vec(), yy.vec());
  }
}
//-----------------------------------------------------------------------------
real PETScMatrix::norm(const Norm type) const
{
  real value = 0.0;
  switch ( type )
  {
  case l1:
    MatNorm(A, NORM_1, &value);
    break;
  case linf:
    MatNorm(A, NORM_INFINITY, &value);
    break;
  case frobenius:
    MatNorm(A, NORM_FROBENIUS, &value);
    break;
  default:
    error("Unknown norm type.");
  }
  return value;
}
//-----------------------------------------------------------------------------
void PETScMatrix::apply(FinalizeType finaltype)
{
  if (finaltype == FINALIZE or finaltype == PETSC_HACK) {
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  } else if ( finaltype == FLUSH ) {
    MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);
  }

}
//-----------------------------------------------------------------------------
void PETScMatrix::zero()
{
  MatZeroEntries(A);
}
//-----------------------------------------------------------------------------
const PETScMatrix& PETScMatrix::operator+= (const PETScMatrix& A)
{
  dolfin_assert(A);
  MatAXPY(this->A, 1.0, A.A, SAME_NONZERO_PATTERN);
  return *this;
}
//-----------------------------------------------------------------------------
const PETScMatrix& PETScMatrix::operator*= (real a)
{
  dolfin_assert(A);
  MatScale(A, a);
  return *this;
}
//-----------------------------------------------------------------------------
const PETScMatrix& PETScMatrix::operator/= (real a)
{
  dolfin_assert(A);
  MatScale(A, 1.0 / a);
  return *this;
}
//-----------------------------------------------------------------------------
const GenericMatrix& PETScMatrix::operator= (const GenericMatrix& A)
{
  //  (*this).down_cast<PETScMatrix>() =  A.down_cast<PETScMatrix>();
  //MatCopy(A.down_cast<PETScMatrix>().A, this->A, DIFFERENT_NONZERO_PATTERN);
  MatCopy(A.down_cast<PETScMatrix>().A, this->A, SAME_NONZERO_PATTERN);
  return *this;
}
//-----------------------------------------------------------------------------
const PETScMatrix& PETScMatrix::operator= (const PETScMatrix& A)
{

  //MatCopy(A.A, (this->A), DIFFERENT_NONZERO_PATTERN);
  MatCopy(A.A, (this->A), SAME_NONZERO_PATTERN);
  return *this;
}
//-----------------------------------------------------------------------------
void PETScMatrix::dup(GenericMatrix& A) 
{
  //  MatDuplicate(A.down_cast<PETScMatrix>().A, MAT_DO_NOT_COPY_VALUES, &this->A);
  MatDuplicate(A.down_cast<PETScMatrix>().A, MAT_COPY_VALUES, &this->A);
}
//-----------------------------------------------------------------------------
PETScMatrix::Type PETScMatrix::type() const
{
  return _type;
}
//-----------------------------------------------------------------------------
void PETScMatrix::disp(uint precision) const
{
  // FIXME: Maybe this could be an option?
  if(MPI::numProcesses() > 1)
    MatView(A, PETSC_VIEWER_STDOUT_WORLD);
  else {
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);   
#ifdef BLOCKED
    PetscObjectSetName((PetscObject) A, " Ab");
#else
    PetscObjectSetName((PetscObject) A, " Anb");
#endif
    MatView(A, PETSC_VIEWER_STDOUT_SELF);
  }


/*
  const uint M = size(0);
  const uint N = size(1);

  // Sparse output
  for (uint i = 0; i < M; i++)
  {
    std::stringstream line;
    line << std::setiosflags(std::ios::scientific);
    line << std::setprecision(precision);
    
    line << "|";
    
    if ( sparse )
    {
      int ncols = 0;
      const int* cols = 0;
      const double* vals = 0;
      MatGetRow(A, i, &ncols, &cols, &vals);
      for (int pos = 0; pos < ncols; pos++)
      {
	       line << " (" << i << ", " << cols[pos] << ", " << vals[pos] << ")";
      }
      MatRestoreRow(A, i, &ncols, &cols, &vals);
    }
    else
    {
      for (uint j = 0; j < N; j++)
      {
        real value = get(i, j);
        if ( fabs(value) < DOLFIN_EPS )
        value = 0.0;	
        line << " " << value;
      }
    }

    line << "|";
    cout << line.str().c_str() << endl;
  }
*/
}
//-----------------------------------------------------------------------------
LinearAlgebraFactory& PETScMatrix::factory() const
{
  return PETScFactory::instance();
}
//-----------------------------------------------------------------------------
Mat PETScMatrix::mat() const
{
  return A;
}
//-----------------------------------------------------------------------------
void PETScMatrix::setType() 
{
  MatType mat_type = getPETScType();
  MatSetType(A, mat_type);
}
//-----------------------------------------------------------------------------
void PETScMatrix::checkType()
{
  switch ( _type )
  {
  case default_matrix:
    return;
  case spooles:
    #if !PETSC_HAVE_SPOOLES
      warning("PETSc has not been complied with Spooles. Using default matrix type");
      _type = default_matrix;
    #endif
    return;
  case superlu:
    #if !PETSC_HAVE_SUPERLU
      warning("PETSc has not been complied with Super LU. Using default matrix type");
      _type = default_matrix;
    #endif
    return;
  case umfpack:
    #if !PETSC_HAVE_UMFPACK
      warning("PETSc has not been complied with UMFPACK. Using default matrix type");
      _type = default_matrix;
    #endif
    return;
  default:
    warning("Requested matrix type unknown. Using default.");
    _type = default_matrix;
  }
}
//-----------------------------------------------------------------------------
MatType PETScMatrix::getPETScType() const
{
  switch ( _type )
  {
  case default_matrix:
    if (MPI::numProcesses() > 1)
      return MATMPIAIJ;
    else
      return MATSEQAIJ;
  case spooles:
      return MATSEQAIJSPOOLES;
  case superlu:
      return MATSUPERLU;
  case umfpack:
      return MATUMFPACK;
  default:
    return "default";
  }
}
//-----------------------------------------------------------------------------
LogStream& dolfin::operator<< (LogStream& stream, const PETScMatrix& A)
{
  stream << "[ PETSc matrix of size " << A.size(0) << " x " << A.size(1) << " ]";
  return stream;
}
//-----------------------------------------------------------------------------

#endif
