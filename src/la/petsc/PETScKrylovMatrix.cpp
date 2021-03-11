// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifdef HAVE_PETSC

#include <dolfin/la/petsc/PETScKrylovMatrix.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/la/petsc/PETScVector.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/MPI.h>

#include <iostream>

using namespace dolfin;

// Mult function
// FIXME: Add an explanation why this function is needed
namespace dolfin
{

  auto usermult(Mat A, Vec x, Vec y) -> int
  {
    void* ctx = nullptr;
    MatShellGetContext(A, &ctx);
    PETScVector xx(x), yy(y);
    ((PETScKrylovMatrix*) ctx)->mult(xx, yy);
    return 0;
  }

}

//-----------------------------------------------------------------------------

PETScKrylovMatrix::PETScKrylovMatrix(const PETScVector& x, const PETScVector& y)
  : A( nullptr )
{
  // Create PETSc matrix
  init(x, y);
}

//-----------------------------------------------------------------------------

PETScKrylovMatrix::~PETScKrylovMatrix()
{
  // Free memory of matrix
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  if ( A ) MatDestroy(&A);
#else
  if ( A ) MatDestroy(A);
#endif
}

//-----------------------------------------------------------------------------

void PETScKrylovMatrix::init(const PETScVector& x, const PETScVector& y)
{
  // Get size and local size of given vector
  int m(0), n(0), M(0), N(0);
  VecGetLocalSize(y.vec(), &m);
  VecGetLocalSize(x.vec(), &n);
  VecGetSize(y.vec(), &M);
  VecGetSize(x.vec(), &N);

  // Free previously allocated memory if necessary
  if ( A )
  {
    // Get size and local size of existing matrix
    int mm(0), nn(0), MM(0), NN(0);
    MatGetLocalSize(A, &mm, &nn);
    MatGetSize(A, &MM, &NN);

    if ( mm == m && nn == n && MM == M && NN == N )
      return;
    else
    {
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
      MatDestroy(&A);
#else
      MatDestroy(A);
#endif
    }
  }

#ifdef DOLFIN_HAVE_MPI
  MatCreateShell(MPI::DOLFIN_COMM, m, n, M, N, (void*) this, &A);
#else
  MatCreateShell(PETSC_COMM_SELF, m, n, M, N, (void*) this, &A);
#endif
  MatShellSetOperation(A, MATOP_MULT, (void (*)()) usermult);
}

//-----------------------------------------------------------------------------

void PETScKrylovMatrix::init(int M, int N)
{
  // Put here to set up arbitrary Shell of global size M,N.
  // Analagous to the matrix being on one processor.

  // Free previously allocated memory if necessary
  if ( A )
    {
      // Get size and local size of existing matrix
      int MM(0), NN(0);
      MatGetSize(A, &MM, &NN);

      if ( MM == M && NN == N )
	return;
      else
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
	MatDestroy(&A);
#else
	MatDestroy(A);
#endif
    }

#ifdef DOLFIN_HAVE_MPI
  MatCreateShell(MPI::DOLFIN_COMM, M, N, M, N, (void*) this, &A);
#else
  MatCreateShell(PETSC_COMM_SELF, M, N, M, N, (void*) this, &A);
#endif
  MatShellSetOperation(A, MATOP_MULT, (void (*)()) usermult);
}

//-----------------------------------------------------------------------------

auto PETScKrylovMatrix::size(size_t dim) const -> dolfin::size_t
{
  int M = 0;
  int N = 0;
  MatGetSize(A, &M, &N);
  dolfin_assert(M >= 0);
  dolfin_assert(N >= 0);

  return (dim == 0 ? static_cast<size_t>(M) : static_cast<size_t>(N));
}

//-----------------------------------------------------------------------------

auto PETScKrylovMatrix::mat() const -> Mat
{
  return A;
}

//-----------------------------------------------------------------------------

void PETScKrylovMatrix::disp(bool, int) const
{
  // Since we don't really have the matrix, we create the matrix by
  // performing multiplication with unit vectors. Used only for debugging.

  warning("Display of PETScKrylovMatrix needs to be fixed.");

/*
  size_t M = size(0);
  size_t N = size(1);
  PETScVector x(N), y(M);
  PETScMatrix A(M, N);


  x = 0.0;
  for (unsigned int j = 0; j < N; j++)
  {
    x(j) = 1.0;
    mult(x, y);
    for (unsigned int i = 0; i < M; i++)
    {
      const real value = y(i);
      if ( fabs(value) > DOLFIN_EPS )
	      A(i, j) = value;
    }
    x(j) = 0.0;
  }

  A.disp(sparse, precision);
*/
}

//-----------------------------------------------------------------------------

auto dolfin::operator<< (LogStream& stream, const PETScKrylovMatrix& A) -> LogStream&
{

#if PETSC_VERSION_MAJOR > 2
#if PETSC_VERSION_MINOR > 3
  MatType type = nullptr;
#else
  const MatType type = nullptr;
#endif
#else
  MatType type = nullptr;
#endif
  MatGetType(A.mat(), &type);
  int m = A.size(0);
  int n = A.size(1);
  stream << "[ PETSc matrix (type " << type << ") of size "
	 << m << " x " << n << " ]";

  return stream;

}

//-----------------------------------------------------------------------------

#endif
