// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//

#include <string>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/Array.h>
#include "JANPACKFactory.h"
#include "JANPACKMat.h"
#include "GenericSparsityPattern.h"

#ifdef HAS_JANPACK

using namespace dolfin;

//-----------------------------------------------------------------------------
JANPACKMat::JANPACKMat():
    Variable("A", "JANPACK matrix"),
    A(0), is_view(false)
{
  // TODO: call JANPACK_Init or something?
}
//-----------------------------------------------------------------------------
JANPACKMat::JANPACKMat(uint M, uint N):
    Variable("A", "JANPACK matrix"),
    A(0), is_view(false)
{
  // TODO: call JANPACK_Init or something?
  // Create JANPACK matrix
  init(M, N);
}
//-----------------------------------------------------------------------------
JANPACKMat::JANPACKMat(const JANPACKMat& A):
  Variable("A", "JANPACK matrix"),
  A(0), is_view(true)
{
  error("Not implemented.");
}
//-----------------------------------------------------------------------------
JANPACKMat::~JANPACKMat()
{
  // Free memory of matrix
  if(A)
    mat_free_crs(A);  

  //  if (!is_view) delete A;
}
//-----------------------------------------------------------------------------
void JANPACKMat::init(uint M, uint N)
{
  // Free previously allocated memory if necessary
  //  if (A) delete A;

  //  message("Calling init %d %d", M, N);
  A = &_A;
  init_mat_crs(A,M,N,50);
  // Not yet implemented
  //  error("JANPACKMat::init(uint, unit) not yet implemented.");
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
  int M = A->M;
  int N = A->N;
  return (dim == 0 ? M : N);
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
  error("Not implemented (set).");
}
//-----------------------------------------------------------------------------
void JANPACKMat::add(const real* block,
		       uint m, const uint* rows,
		       uint n, const uint* cols)
{
  dolfin_assert(A); 
  const real *bp = &block[0];
  for(uint i = 0 ; i < m; i++)
    for(uint j = 0; j < n; j++)
      mat_add_crs(A, rows[i], cols[j], *(bp++));
  
  //error("Not implemented. (add)");
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
  mat_zero_crs(A);
  //  error("Not implemented. (zero)");
}
//-----------------------------------------------------------------------------
void JANPACKMat::apply(FinalizeType finaltype)
{
  finalize_mat_crs(A);
  //  error("Not implemented. (apply)");
}
//-----------------------------------------------------------------------------
void JANPACKMat::disp(uint precision) const
{
  dolfin_assert(A); 
  error("Not implemented.");
}
//-----------------------------------------------------------------------------
void JANPACKMat::ident(uint m, const uint* rows)
{
  dolfin_assert(A); 

  for(uint i = 0; i < m; i ++) {
    mat_zero_row_crs(A, rows[i]);
    mat_insert_crs(A, rows[i], rows[i], 1.0);
  }

  //  error("Not implemented.");
}
//-----------------------------------------------------------------------------
void JANPACKMat::zero(uint m, const uint* rows)
{
  dolfin_assert(A); 
  error("Not implemented.");
}
//-----------------------------------------------------------------------------
void JANPACKMat::mult(const GenericVector& x_, GenericVector& Ax_, bool transposed) const
{
  dolfin_assert(A); 
  error("Not implemented.");
}
//-----------------------------------------------------------------------------
void JANPACKMat::getrow(uint row, Array<uint>& columns, Array<real>& values) const
{
  dolfin_assert(A); 
  error("Not implemented.");
}
//-----------------------------------------------------------------------------
void JANPACKMat::setrow(uint row, const Array<uint>& columns, const Array<real>& values)
{
  error("Not implemented.");
}
//-----------------------------------------------------------------------------
LinearAlgebraFactory& JANPACKMat::factory() const
{
  return JANPACKFactory::instance();
}
//-----------------------------------------------------------------------------
Mat_crs *JANPACKMat::mat() const
{
  return A;
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


#endif
