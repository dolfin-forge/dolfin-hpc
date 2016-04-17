// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-05-15
// Last changed: 2008-05-15

#include <dolfin/common/Array.h>
#include <dolfin/la/LinearAlgebraFactory.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/SingularSolver.h>
#include <dolfin/main/MPI.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
SingularSolver::SingularSolver(SolverType solver_type,
                               PreconditionerType pc_type) :
    Parametrized(),
    linear_solver(solver_type, pc_type),
    B(0),
    y(0),
    c(0)
{
  // Set parameters for linear solver
  linear_solver.set("parent", *this);

  //linear_solver.set("Krylov monitor convergence", true);
}
//-----------------------------------------------------------------------------
SingularSolver::~SingularSolver()
{
  delete B;
  delete y;
  delete c;
}
//-----------------------------------------------------------------------------
dolfin::uint SingularSolver::solve(const GenericMatrix& A, GenericVector& x,
                                   const GenericVector& b)
{
  message("Solving singular system...");

  // Initialize data structures for extended system
  init(A);

  // Create extended system
  create(A, b, 0);

  // Solve extended system
  const uint num_iterations = linear_solver.solve(*B, *y, *c);

  // Extract solution
  x.init(y->size() - 1);
  real* vals = new real[y->size()];
  y->get(vals);
  x.set(vals);
  delete[] vals;

  return num_iterations;
}
//-----------------------------------------------------------------------------
dolfin::uint SingularSolver::solve(const GenericMatrix& A, GenericVector& x,
                                   const GenericVector& b,
                                   const GenericMatrix& M)
{
  message("Solving singular system...");

  // Initialize data structures for extended system
  init(A);

  // Create extended system
  create(A, b, &M);

  // Solve extended system
  const uint num_iterations = linear_solver.solve(*B, *y, *c);

  // Extract solution
  x.init(y->size() - 1);
  real* vals = new real[y->size()];
  y->get(vals);
  x.set(vals);
  delete[] vals;

  return num_iterations;
}
//-----------------------------------------------------------------------------
void SingularSolver::init(const GenericMatrix& A)
{
  // Check size of system
  if (A.size(0) != A.size(1)) error("Matrix must be square.");
  if (A.size(0) == 0) error("Matrix size must be non-zero.");

  // Get dimension
  const uint N = A.size(0);

  // Check if we have already initialized system
  if (B && B->size(0) == N + 1 && B->size(1) == N + 1) return;

  // Delete any old data
  delete B;
  delete y;
  delete c;

  // Create sparsity pattern for B
  uint dims[2] = { N + 1, N + 1 };
  if (MPI::numProcesses() > 1)
  {
    error("SingularSolver : not implemented in parallel");
  };
  SparsityPattern s(2, dims);

  // Copy sparsity pattern for A and last column
  Array<uint> columns;
  Array<real> dummy;
  uint num_rows[2] = { 1, 0 };
  uint * rows[2];
  rows[0] = new uint[1];
  rows[0][0] = 0;
  rows[1] = new uint[A.size(0) + 1];
  std::fill_n(rows[1], A.size(0) + 1, 0.0);
  for (rows[0][0] = 0; rows[0][0] < N; ++rows[0][0])
  {
    // Get row
    A.getrow(rows[0][0], columns, dummy);

    // Copy columns to array
    num_rows[1] = columns.size() + 1;
    for (uint j = 0; j < columns.size(); ++j)
    {
      rows[1][j] = columns[j];
    }

    // Add last entry
    rows[1][num_rows[1] - 1] = N;

    // Insert into sparsity pattern
    s.insert(num_rows, rows);
  }

  // Add last row
  rows[0][0] = N;
  num_rows[1] = N;
  for (uint j = 0; j < num_rows[1]; ++j)
  {
    rows[1][j] = j;
  }
  // Insert into sparsity pattern
  s.insert(num_rows, rows);

  delete rows[0];
  delete rows[1];

  // Create matrix and vector
  B = A.factory().createMatrix();
  y = A.factory().createVector();
  c = A.factory().createVector();
  B->init(s);
  y->init(N + 1);
  c->init(N + 1);
}
//-----------------------------------------------------------------------------
void SingularSolver::create(const GenericMatrix& A, const GenericVector& b,
                            const GenericMatrix* M)
{
  dolfin_assert(B);
  dolfin_assert(c);

  message("Creating extended hopefully non-singular system...");

  // Reset matrix
  B->zero();

  // Copy rows from A into B
  const uint N = A.size(0);
  Array<uint> columns;
  Array<real> values;
  for (uint i = 0; i < N; i++)
  {
    A.getrow(i, columns, values);
    B->setrow(i, columns, values);
  }

  // Compute lumped mass matrix
  columns.resize(N);
  values.resize(N);
  if (M)
  {
    GenericVector* ones = A.factory().createVector();
    GenericVector* z = A.factory().createVector();
    ones->init(N);
    z->init(N);
    *ones = 1.0;
    M->mult(*ones, *z);
    for (uint i = 0; i < N; i++)
    {
      columns[i] = i;
      values[i] = (*z)[i];
    }
    delete ones;
    delete z;
  }
  else
  {
    for (uint i = 0; i < N; i++)
    {
      columns[i] = i;
      values[i] = 1.0;
    }
  }

  // Add last row
  B->setrow(N, columns, values);

  // Add last column
  for (uint i = 0; i < N; i++)
    B->set(&values[i], 1, &i, 1, &N);

  // Copy values from b into c
  real* vals = new real[N + 1];
  b.get(vals);
  vals[N] = 0.0;
  c->set(vals);
  delete[] vals;

  // Apply changes
  B->apply();
  c->apply();
}
//-----------------------------------------------------------------------------
