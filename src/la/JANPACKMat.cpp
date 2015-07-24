// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//

#include <string>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/Array.h>
#include <dolfin/la/JANPACKFactory.h>
#include <dolfin/la/JANPACKMat.h>
#include <dolfin/la/JANPACKVec.h>
#include <dolfin/la/GenericSparsityPattern.h>


#ifdef HAVE_JANPACK


#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_MPI
#include <dolfin/main/MPI.h>
#endif

#include <janpack/spmv.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
JANPACKMat::JANPACKMat():
    Variable("A", "JANPACK matrix"), is_view(false)
#ifdef HAVE_JANPACK_MPI
    , A(&AA)
#endif
{
  // TODO: call JANPACK_Init or something?
}
//-----------------------------------------------------------------------------
JANPACKMat::JANPACKMat(uint M, uint N):
    Variable("A", "JANPACK matrix"),
    is_view(false)
 {
  // TODO: call JANPACK_Init or something?
  // Create JANPACK matrix
  init(M, N);
}
//-----------------------------------------------------------------------------
JANPACKMat::JANPACKMat(const JANPACKMat& A):
  Variable("A", "JANPACK matrix"),
   is_view(true)
{
  error("Not implemented");
}
//-----------------------------------------------------------------------------
JANPACKMat::~JANPACKMat()
{
  // Free memory of matrix
  //  if(A)
    jp_mat_free(A);

  //  if (!is_view) delete A;
}
//-----------------------------------------------------------------------------
void JANPACKMat::init(uint M, uint N)
{
  // Free previously allocated memory if necessary
  //  if (A) delete A;

  //  A = &AA;
  jp_mat_init(A, M, N);

  //  jp_mat_setopt(A, JP_MAT_SORTED);
  // Not yet implemented
  //  error("JANPACKMat::init(uint, unit) not yet implemented.");
}
//-----------------------------------------------------------------------------
void JANPACKMat::init(uint M, uint N, bool distributed)
{
  init(M, N);
}
//-----------------------------------------------------------------------------
void JANPACKMat::init(const GenericSparsityPattern& sparsity_pattern)
{

  const SparsityPattern& spattern =
    reinterpret_cast<const SparsityPattern&>(sparsity_pattern);

  init(spattern.size(0), spattern.size(1));

  // error("Not implemented.  (init)");
}
//-----------------------------------------------------------------------------
JANPACKMat* JANPACKMat::copy() const
{

  error("JANPACKMat::copy not yet implemented.");

  JANPACKMat *mcopy = new JANPACKMat();

  return mcopy;
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKMat::size(uint dim) const
{
  dolfin_assert(A);
  uint32_t M = 0;
  uint32_t N = 0;
  jp_mat_size(const_cast<jp_mat_type *>(A), &M, &N);
  return (dim == 0 ? M : N);
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKMat::nz() const
{
  return jp_mat_nz(const_cast<jp_mat_type *>(A), 1);
}
//-----------------------------------------------------------------------------
void JANPACKMat::get(real* block,
		       uint m, const uint* rows,
		       uint n, const uint* cols) const
{
  dolfin_assert(A);
  // for each row in rows
  //A->ExtractGlobalRowCopy(...)

  // Not yet implemented
  error("JANPACKMat::get not yet implemented.");
}
//-----------------------------------------------------------------------------
void JANPACKMat::set(const real* block,
		       uint m, const uint* rows,
		       uint n, const uint* cols)
{
  dolfin_assert(A);

  jp_mat_set_block(A,
	       m, const_cast<uint*>(rows),
	       n, const_cast<uint*>(cols),
	       const_cast<real*>(block));
}
//-----------------------------------------------------------------------------
void JANPACKMat::add(const real* block,
		       uint m, const uint* rows,
		       uint n, const uint* cols)
{
  dolfin_assert(A);

  jp_mat_add_block(A,
	       m, const_cast<uint*>(rows),
	       n, const_cast<uint*>(cols),
	       const_cast<real*>(block));
}
//-----------------------------------------------------------------------------
real JANPACKMat::norm(std::string norm_type) const
{
  error("Not implemented.");
  return 0.0;
}
//-----------------------------------------------------------------------------
void JANPACKMat::zero()
{
  dolfin_assert(A);
  jp_mat_zero(A);
}
//-----------------------------------------------------------------------------
void JANPACKMat::apply(FinalizeType finaltype)
{

  jp_mat_finalize(A);

}
//-----------------------------------------------------------------------------
void JANPACKMat::disp(uint precision) const
{
#ifndef HAVE_JANPACK_MPI
  jp_mat_print(const_cast<jp_mat_type *>(A));
#endif
}
//-----------------------------------------------------------------------------
void JANPACKMat::ident(uint m, const uint* rows)
{
  dolfin_assert(A);

  for(uint i = 0; i < m; i ++) {
    //jp_mat_ident(A, rows[i]);
    jp_mat_zero_row(A, rows[i]);
    jp_mat_insert(A, rows[i], rows[i], 1.0);
  }
}
//-----------------------------------------------------------------------------
void JANPACKMat::zero(uint m, const uint* rows)
{
  dolfin_assert(A);

  for (uint i = 0; i < m; i++)
    jp_mat_zero_row(A,rows[i]);
}
//-----------------------------------------------------------------------------
void JANPACKMat::mult(const GenericVector& x, GenericVector& y, bool transposed) const
{
  dolfin_assert(A);
  const JANPACKVec& xx = x.down_cast<JANPACKVec>();
  JANPACKVec& yy = y.down_cast<JANPACKVec>();
  if (transposed)
    yy.init(xx.local_size());
  else
    yy.init(xx.local_size());

  jp_spmv(const_cast<jp_mat_type *>(A), xx.vec(), yy.vec());
}
//-----------------------------------------------------------------------------
void JANPACKMat::getrow(uint row, Array<uint>& columns, Array<real>& values) const
{

  uint  n, *c = 0;
  real *v = 0;

  c = new uint[size(0)];
  v = new real[size(0)];

  jp_mat_getrow(const_cast<jp_mat_type *>(A), row, c, v, &n);

  for (uint i = 0; i < n; i++)
  {
    columns.push_back(c[i]);
    values.push_back(v[i]);
  }

  delete c;
  delete v;

}
//-----------------------------------------------------------------------------
void JANPACKMat::setrow(uint row, const Array<uint>& columns,
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
LinearAlgebraFactory& JANPACKMat::factory() const
{
  return JANPACKFactory::instance();
}
//-----------------------------------------------------------------------------
jp_mat_type* JANPACKMat::mat() const
{
  return const_cast<jp_mat_type *>(A);
}
//-----------------------------------------------------------------------------
const JANPACKMat& JANPACKMat::operator*= (real a)
{
  dolfin_assert(A);
  error("Not implemented.");
  return *this;
}
//-----------------------------------------------------------------------------
const JANPACKMat& JANPACKMat::operator/= (real a)
{
  dolfin_assert(A);
  error("Not implemented.");
  return *this;
}
//-----------------------------------------------------------------------------
const GenericMatrix& JANPACKMat::operator= (const GenericMatrix& A)
{
  //  jp_mat_copy(A.down_cast<JANPACKMat>().A, this->A);
  error("Not implemented.");
  return *this;
}
//-----------------------------------------------------------------------------
void JANPACKMat::dup(const JANPACKMat& A)
{
  uint range[2];
  uint m;

  jp_mat_range(const_cast<jp_mat_type *>(A.mat()), &range[0], &range[1]);
  m = range[1] - range[0];
  init(m, m);
  //  error("Not implemented.");
  //  jp_mat_dup_crs(, this->A);
}
//-----------------------------------------------------------------------------
#endif
