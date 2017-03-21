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

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include <dolfin/common/Array.h>
#include <dolfin/common/system.h>
#include <dolfin/la/PETScVector.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/GenericSparsityPattern.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/PETScFactory.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

#include <iostream>
#include <sstream>
#include <iomanip>

namespace dolfin
{

//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix() :
    Variable("A", "a sparse matrix"),
    A(NULL),
    AA_sub(NULL),
    is_distributed_(false),
    rstart_(0),
    rend_(0)
{
}
//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix(Mat A) :
    Variable("A", "a sparse matrix"),
    A(A),
    AA_sub(NULL),
    is_distributed_(false),
    rstart_(0),
    rend_(0)
{
}
//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix(uint M, uint N, bool distributed) :
    Variable("A", "a sparse matrix"),
    A(NULL),
    AA_sub(NULL),
    is_distributed_(false),
    rstart_(0),
    rend_(0)
{
  init(M, N, distributed);
}
//-----------------------------------------------------------------------------
PETScMatrix::PETScMatrix(const PETScMatrix& A) :
    Variable("A", "PETSc matrix"),
    A(NULL),
    AA_sub(NULL),
    is_distributed_(false),
    rstart_(0),
    rend_(0)
{
  *this = A;
}
//-----------------------------------------------------------------------------
PETScMatrix::~PETScMatrix()
{
  clear();
}
//-----------------------------------------------------------------------------
void PETScMatrix::clear()
{
  if (A)
  {
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    MatDestroy(&A);
#else
    MatDestroy(A);
#endif
  }

  if (AA_sub)
  {
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    MatDestroyMatrices(1, &AA_sub);
#else
    MatDestroyMatrices(1, AA_sub);
#endif
    AA_sub = NULL;
  }
}
//-----------------------------------------------------------------------------
void PETScMatrix::init(uint M, uint N)
{
  init(M, N, true);
}
//-----------------------------------------------------------------------------
void PETScMatrix::init(uint M, uint N, bool distributed)
{
  // Free previously allocated memory if necessary
  clear();

  // Create a sparse matrix in compressed row format
  if (dolfin::MPI::size() > 1 && distributed)
  {
    is_distributed_ = true;
    // Create PETSc parallel matrix with a guess for number of diagonal and
    /// number of off-diagonal non-zeroes.
    // Note that guessing too high leads to excessive memory usage.
#ifdef HAVE_MPI
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3
    MatCreateAIJ(dolfin::MPI::DOLFIN_COMM, M, N, PETSC_DETERMINE,
                 PETSC_DETERMINE, 120, PETSC_NULL, 120, PETSC_NULL, &A);
#else
    MatCreateMPIAIJ(dolfin::MPI::DOLFIN_COMM, M, N, PETSC_DETERMINE,
                    PETSC_DETERMINE, 120, PETSC_NULL, 120, PETSC_NULL, &A);
#endif
#endif
  }
  else
  {
    // Create PETSc sequential matrix with a guess for number of non-zeroes
    MatCreateSeqAIJ(PETSC_COMM_SELF, M, N, 50, PETSC_NULL, &A);

#if PETSC_VERSION_MAJOR > 2
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 0
    MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
#else
    MatSetOption(A, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE);
#endif
#else
    MatSetOption(A, MAT_KEEP_ZEROED_ROWS);
#endif
    MatSetFromOptions(A);
  }

  MatGetOwnershipRange(A, &rstart_, &rend_);
}
//-----------------------------------------------------------------------------
void PETScMatrix::init(uint M, uint N, uint const* nz)
{
  dolfin_assert(is_distributed_ == false);

  // Free previously allocated memory if necessary
  clear();

  // Create a sparse matrix in compressed row format
  MatCreate(PETSC_COMM_SELF, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, M, N);
  MatSetType(A, MATSEQAIJ);

#if PETSC_VERSION_MAJOR > 2
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 0
  MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
#else
  MatSetOption(A, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE);
#endif
#else
  MatSetOption(A, MAT_KEEP_ZEROED_ROWS);
#endif
  MatSetFromOptions(A);
  MatSeqAIJSetPreallocation(A, PETSC_DEFAULT, (int*) nz);
  MatZeroEntries(A);

  MatGetOwnershipRange(A, &rstart_, &rend_);
}
//-----------------------------------------------------------------------------
void PETScMatrix::init(uint M, uint N, uint const* d_nzrow, uint const* o_nzrow)
{
  dolfin_assert(is_distributed_ == true);

  // Free previously allocated memory if necessary
  clear();

  // Create PETSc parallel matrix with a guess for number of diagonal (50 in this case)
  // and number of off-diagonal non-zeroes (50 in this case).
  // Note that guessing too high leads to excessive memory usage.
  // In order to not waste any memory one would need to specify d_nnz and o_nnz.
#ifdef HAVE_MPI
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3
  MatCreateAIJ(MPI::DOLFIN_COMM, M, N, PETSC_DETERMINE, PETSC_DETERMINE,
               PETSC_DETERMINE, (PetscInt*) d_nzrow, PETSC_DETERMINE, (PetscInt*) o_nzrow,
               &A);
#else
  MatCreateMPIAIJ(MPI::DOLFIN_COMM, M, N, PETSC_DETERMINE, PETSC_DETERMINE,
                  PETSC_DETERMINE, (int*)d_nzrow, PETSC_DETERMINE,
                  (PetscInt*)o_nzrow, &A);
#endif
#endif

#if PETSC_VERSION_MAJOR > 2
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 0
  MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
#else
  MatSetOption(A, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE);
#endif
#else
  MatSetOption(A, MAT_KEEP_ZEROED_ROWS);
#endif
  MatSetFromOptions(A);
  MatZeroEntries(A);

  MatGetOwnershipRange(A, &rstart_, &rend_);
}
//-----------------------------------------------------------------------------
void PETScMatrix::init(const GenericSparsityPattern& sparsity_pattern)
{
  if (sparsity_pattern.is_distributed())
  {
    is_distributed_ = true;
    uint p = dolfin::MPI::rank();
    SparsityPattern const& spattern =
        reinterpret_cast<SparsityPattern const&>(sparsity_pattern);
    uint local_size = spattern.range_size(p);
    uint* d_nzrow = new uint[local_size];
    uint* o_nzrow = new uint[local_size];
    spattern.numNonZeroPerRow(p, d_nzrow, o_nzrow);
    uint spattern_size_0 = spattern.size(0);
    uint spattern_size_1 = spattern.size(1);
    const_cast<SparsityPattern&>(spattern).clear();
    init(spattern_size_0, spattern_size_1, d_nzrow, o_nzrow);
    delete[] d_nzrow;
    delete[] o_nzrow;
  }
  else
  {
    SparsityPattern const& spattern =
        reinterpret_cast<SparsityPattern const&>(sparsity_pattern);
    uint* nzrow = new uint[spattern.size(0)];
    spattern.numNonZeroPerRow(nzrow);
    init(spattern.size(0), spattern.size(1), nzrow);
    delete[] nzrow;
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
dolfin::uint PETScMatrix::nz() const
{
  MatInfo info;
  MatGetInfo(A, MAT_GLOBAL_SUM, &info);

  return info.nz_used;
}
//-----------------------------------------------------------------------------
void PETScMatrix::get(real* block, uint m, uint const* rows, uint n,
                      uint const* cols) const
{
  dolfin_assert(A);
  MatGetValues(A, static_cast<int>(m),
               reinterpret_cast<int*>(const_cast<uint*>(rows)),
               static_cast<int>(n),
               reinterpret_cast<int*>(const_cast<uint*>(cols)), block);
}
//-----------------------------------------------------------------------------
void PETScMatrix::set(real const* block, uint m, uint const* rows, uint n,
                      uint const* cols)
{
  dolfin_assert(A);
  MatSetValues(A, static_cast<int>(m),
               reinterpret_cast<int*>(const_cast<uint*>(rows)),
               static_cast<int>(n),
               reinterpret_cast<int*>(const_cast<uint*>(cols)), block,
               INSERT_VALUES);
}
//-----------------------------------------------------------------------------
void PETScMatrix::add(real const* block, uint m, uint const* rows, uint n,
                      uint const* cols)
{
  dolfin_assert(A);
  MatSetValues(A, static_cast<int>(m),
               reinterpret_cast<int*>(const_cast<uint*>(rows)),
               static_cast<int>(n),
               reinterpret_cast<int*>(const_cast<uint*>(cols)), block,
               ADD_VALUES);
}
//-----------------------------------------------------------------------------
real PETScMatrix::norm(std::string norm_type) const
{
  real norm;
  if (norm_type == "l1")
  {
    MatNorm(A, NORM_1, &norm);
  }
  else if (norm_type == "linf")
  {
    MatNorm(A, NORM_INFINITY, &norm);
  }
  else if (norm_type == "frobenius")
  {
    MatNorm(A, NORM_FROBENIUS, &norm);
  }
  else
  {
    error("Unknown norm type in uPETScMatrix.");
    return 0.0;
  }
  return norm;
}
//-----------------------------------------------------------------------------
void PETScMatrix::getrow(uint row, Array<uint>& columns,
                         Array<real>& values) const
{
  const PetscInt *cols = NULL;
  const PetscScalar *vals = NULL;
  PetscInt ncols = 0;
  if (row >= static_cast<uint>(rstart_) && row < static_cast<uint>(rend_))
  {
    // Local row
    MatGetRow(A, row, &ncols, &cols, &vals);
    columns.assign(reinterpret_cast<uint*>(const_cast<int*>(cols)),
                   reinterpret_cast<uint*>(const_cast<int*>(cols + ncols)));
    values.assign(vals, vals + ncols);
    MatRestoreRow(A, row, &ncols, &cols, &vals);
  }
  else
  {
    // Non-local row
    if (AA_sub)
    {
      error("Trying to get ghosted row %d but no SubMatrices defined.\n"
            "Ownership range %d, %d",
            row, rstart_, rend_);
    }

    std::map<int, int>::const_iterator it = mapping_.find(row);
    MatGetRow(AA_sub[0], it->second, &ncols, &cols, &vals);
    columns.assign(reinterpret_cast<uint*>(const_cast<int*>(cols)),
                   reinterpret_cast<uint*>(const_cast<int*>(cols + ncols)));
    values.assign(vals, vals + ncols);
    MatRestoreRow(AA_sub[0], it->second, &ncols, &cols, &vals);
  }
}
//-----------------------------------------------------------------------------
void PETScMatrix::getrows_offproc(std::set<uint> const& rows)
{
  if (!is_distributed_)
  {
    return;
  }

  int *_cols = new int[size(0)];
  int *_rows = new int[rows.size()];

  mapping_.clear();

  uint i = 0;
  for (std::set<uint>::const_iterator it = rows.begin(); it != rows.end(); ++it)
  {
    if (*it < static_cast<uint>(rstart_) && *it >= static_cast<uint>(rend_))
    {
      _rows[i] = *it;
      mapping_[*it] = i++;
    }
  }

  for (uint j = 0; j < size(0); ++j)
  {
    _cols[j] = j;
  }

  IS irow, icol;
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  ISCreateGeneral(PETSC_COMM_SELF, i, _rows, PETSC_COPY_VALUES, &irow);
  ISCreateGeneral(PETSC_COMM_SELF, size(0), _cols, PETSC_COPY_VALUES, &icol);
#else
  ISCreateGeneral(PETSC_COMM_SELF, i, _rows, &irow);
  ISCreateGeneral(PETSC_COMM_SELF, size(0), _cols, &icol);
#endif

  if (!AA_sub)
  {
    MatGetSubMatrices(A, 1, &irow, &icol, MAT_INITIAL_MATRIX, &AA_sub);
  }
  else
  {
    MatGetSubMatrices(A, 1, &irow, &icol, MAT_REUSE_MATRIX, &AA_sub);
  }

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  ISDestroy(&irow);
  ISDestroy(&icol);
#else
  ISDestroy(irow);
  ISDestroy(icol);
#endif

  delete[] _cols;
  delete[] _rows;

}
//-----------------------------------------------------------------------------
void PETScMatrix::setrow(uint row, const Array<uint>& columns,
                         const Array<real>& values)
{
  set(&values[0], 1, &row, columns.size(), &columns[0]);
}
//-----------------------------------------------------------------------------
void PETScMatrix::zero(uint m, uint const* rows)
{
  IS is = 0;
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  ISCreateGeneral(PETSC_COMM_SELF, static_cast<int>(m),
                  reinterpret_cast<int*>(const_cast<uint*>(rows)),
                  PETSC_COPY_VALUES, &is);
#else
  ISCreateGeneral(PETSC_COMM_SELF, static_cast<int>(m),
      reinterpret_cast<int*>(const_cast<uint*>(rows)), &is);
#endif
  PetscScalar null = 0.0;
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  MatZeroRowsIS(A, is, null, PETSC_NULL, PETSC_NULL);
  ISDestroy(&is);
#else
  MatZeroRowsIS(A, is, null);
  ISDestroy(is);
#endif
}
//-----------------------------------------------------------------------------
void PETScMatrix::ident(uint m, uint const* rows)
{
  IS is = 0;
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  ISCreateGeneral(PETSC_COMM_SELF, static_cast<int>(m),
                  reinterpret_cast<int*>(const_cast<uint*>(rows)),
                  PETSC_COPY_VALUES, &is);
#else
  ISCreateGeneral(PETSC_COMM_SELF, static_cast<int>(m),
      reinterpret_cast<int*>(const_cast<uint*>(rows)), &is);
#endif
  PetscScalar one = 1.0;
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  MatZeroRowsIS(A, is, one, PETSC_NULL, PETSC_NULL);
  ISDestroy(&is);
#else
  MatZeroRowsIS(A, is, one);
  ISDestroy(is);
#endif
}
//-----------------------------------------------------------------------------
void PETScMatrix::mult(const GenericVector& x, GenericVector& y,
                       bool transposed) const
{
  const PETScVector& xx = x.down_cast<PETScVector>();
  PETScVector& yy = y.down_cast<PETScVector>();

  if (transposed)
  {
    if (size(0) != xx.size())
    {
      error("Matrix and vector dimensions mismatch for matrix-vector product.");
    }
    yy.init(xx.local_size());
    MatMultTranspose(A, xx.vec(), yy.vec());
  }
  else
  {
    if (size(1) != xx.size())
    {
      error("Matrix and vector dimensions mismatch for matrix-vector product.");
    }
    yy.init(xx.local_size());
    MatMult(A, xx.vec(), yy.vec());
  }
}
//-----------------------------------------------------------------------------
real PETScMatrix::norm(const Norm type) const
{
  real value = 0.0;
  switch (type)
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
      break;
    }
  return value;
}
//-----------------------------------------------------------------------------
void PETScMatrix::apply(FinalizeType finaltype)
{
  if (finaltype == FINALIZE)
  {
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  }
  else if (finaltype == FLUSH)
  {
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
const PETScMatrix& PETScMatrix::operator+=(const PETScMatrix& A)
{
  dolfin_assert(this->A);
  MatAXPY(this->A, 1.0, A.A, SAME_NONZERO_PATTERN);
  return *this;
}
//-----------------------------------------------------------------------------
const PETScMatrix& PETScMatrix::operator*=(real a)
{
  dolfin_assert(A);
  MatScale(A, a);
  return *this;
}
//-----------------------------------------------------------------------------
const PETScMatrix& PETScMatrix::operator/=(real a)
{
  dolfin_assert(A);
  MatScale(A, 1.0 / a);
  return *this;
}
//-----------------------------------------------------------------------------
const GenericMatrix& PETScMatrix::operator=(const GenericMatrix& A)
{
  if (&A != this)
  {
    MatCopy(A.down_cast<PETScMatrix>().A, this->A, SAME_NONZERO_PATTERN);
  }
  return *this;
}
//-----------------------------------------------------------------------------
const PETScMatrix& PETScMatrix::operator=(const PETScMatrix& A)
{
  if (&A != this)
  {
    MatCopy(A.A, (this->A), SAME_NONZERO_PATTERN);
  }
  return *this;
}
//-----------------------------------------------------------------------------
void PETScMatrix::dup(GenericMatrix& A)
{
  MatDuplicate(A.down_cast<PETScMatrix>().A, MAT_COPY_VALUES, &this->A);
}
//-----------------------------------------------------------------------------
void PETScMatrix::disp(uint precision) const
{
  section("PETScMatrix");
  MatInfo info;

  section("Local");
  MatGetInfo(A, MAT_LOCAL, &info);
  print(info);
  endblock();
  skip();

  if (is_distributed_)
  {
    section("Global");
    MatGetInfo(A, MAT_GLOBAL_SUM, &info);
    print(info);
    endblock();
    skip();
  }

  if (logm.getDebugLevel() == 0)
  {
    return;
  }

  if (is_distributed_)
  {

    MatView(A, PETSC_VIEWER_STDOUT_WORLD);
  }
  else
  {
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
#ifdef BLOCKED
    PetscObjectSetName((PetscObject) A, " Ab");
#else
    PetscObjectSetName((PetscObject) A, " Anb");
#endif
    MatView(A, PETSC_VIEWER_STDOUT_SELF);
  }
  endblock();
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
void PETScMatrix::print(MatInfo const& info) const
{
  message("%24s : %lu", "block_size", uidx(info.block_size));
  message("%24s : %lu", "nz_allocated", uidx(info.nz_allocated));
  message("%24s : %lu", "nz_used", uidx(info.nz_used));
  message("%24s : %lu", "nz_unneeded", uidx(info.nz_unneeded));
  message("%24s : %s" , "memory", human_readable(info.memory).c_str());
  message("%24s : %lu", "assemblies", uidx(info.assemblies));
  message("%24s : %lu", "mallocs", uidx(info.mallocs));
  message("%24s : %lu", "fill_ratio_given", uidx(info.fill_ratio_given));
  message("%24s : %lu", "fill_ratio_needed", uidx(info.fill_ratio_needed));
  message("%24s : %lu", "factor_mallocs", uidx(info.factor_mallocs));
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* HAVE_PETSC */
