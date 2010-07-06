// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//

#ifdef HAS_JANPACK

#include <dolfin/math/dolfin_math.h>
#include <dolfin/log/dolfin_log.h>
#include "JANPACKVec.h"
#include "JANPACKFactory.h"
#include <dolfin/main/MPI.h>

#include <dolfin/common/Array.h>

#include <set>
#include <map>

using namespace dolfin;

//-----------------------------------------------------------------------------
JANPACKVec::JANPACKVec():
    Variable("x", "a sparse vector"),
    x(0), is_view(false), is_ghosted(false)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
JANPACKVec::JANPACKVec(uint N):
    Variable("x", "a sparse vector"), 
    x(0), is_view(false), is_ghosted(false)
{
  // Create PETSc vector
  init(N);
}
//-----------------------------------------------------------------------------
JANPACKVec::JANPACKVec(const JANPACKVec& v):
    Variable("x", "a vector"),
    x(0), is_view(false), is_ghosted(false)
{
  *this = v;
}
//-----------------------------------------------------------------------------
JANPACKVec::~JANPACKVec()
{
  if (x && !is_view)
    vec_free(x);
}
//-----------------------------------------------------------------------------
void JANPACKVec::init(uint N)
{
  
  // Two cases:
  //
  //   1. Already allocated and dimension changes -> reallocate
  //   2. Not allocated -> allocate
  //
  // Otherwise do nothing
  
  if ((x && this->size() == N ) || (this->local_size() == N))
  {
    vec_zero(x);
    return;      
  }
  else
  {
    if (x && !is_view) 
    {
      vec_free(x);
      x = 0;
    }
    
  }
  
  // Create vector
  init_vec(&_x, N);
  x = &_x;
}
//-----------------------------------------------------------------------------
JANPACKVec* JANPACKVec::copy() const
{
  JANPACKVec* v = new JANPACKVec(*this); 
  return v; 
}
//-----------------------------------------------------------------------------
void JANPACKVec::get(real* values) const
{
  memcpy(values, x->x, x->n * sizeof(real));
}
//-----------------------------------------------------------------------------
void JANPACKVec::set(real* values)
{
  memcpy(x->x,values, x->n * sizeof(real));
}
//-----------------------------------------------------------------------------
void JANPACKVec::add(real* values)
{
  dolfin_assert(x);
  
  error("Not implemented.");

  /*  int m = static_cast<int>(size());
  int* rows = new int[m];
  for (int i = 0; i < m; i++)
    rows[i] = i;

  VecSetValues(x, m, rows, values, ADD_VALUES);

  delete [] rows;
  */
}
//-----------------------------------------------------------------------------
void JANPACKVec::get(real* block, uint m, const uint* rows) const
{
  dolfin_assert(x);
  for(uint i = 0; i < m ; i++)
    block[i] = x->x[rows[i]];
}
//-----------------------------------------------------------------------------
void JANPACKVec::set(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x);
  for(uint i = 0; i < m ; i++)
    x->x[rows[i]] = block[i];
  //  error("Not implemented.");
}
//-----------------------------------------------------------------------------
void JANPACKVec::add(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x);

  for(uint i = 0; i < m ; i++)
    x->x[rows[i]] += block[i];
  
  //  error("Not implemented.");
}
//-----------------------------------------------------------------------------
void JANPACKVec::apply(FinalizeType finaltype)
{
  
}
//-----------------------------------------------------------------------------
void JANPACKVec::zero()
{
  dolfin_assert(x);
  vec_zero(x);
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKVec::size() const
{
  int n = 0;
  if (x)
    n  = x->m;

  return static_cast<uint>(n);
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKVec::local_size() const
{
  int n = 0;
  if (x) 
    n = x->n;

  return static_cast<uint>(n);
}
//-----------------------------------------------------------------------------
const GenericVector& JANPACKVec::operator= (const GenericVector& v)
{
  *this = v.down_cast<JANPACKVec>();
  return *this; 
}
//-----------------------------------------------------------------------------
const JANPACKVec& JANPACKVec::operator= (const JANPACKVec& v)
{
  dolfin_assert(v.x);

  init(v.local_size());
  //  vec_copy(x, v.x);
  vec_copy(v.x, x);

  return *this; 
}
//-----------------------------------------------------------------------------
const JANPACKVec& JANPACKVec::operator= (real a)
{
  dolfin_assert(x);
  // VecSet(x, a);
  return *this; 
}
//-----------------------------------------------------------------------------
const JANPACKVec& JANPACKVec::operator+= (const GenericVector& x)
{
  this->axpy(1.0, x); 
  return *this;
}
//-----------------------------------------------------------------------------
const JANPACKVec& JANPACKVec::operator-= (const GenericVector& x)
{
  this->axpy(-1.0, x); 
  return *this;
}
//-----------------------------------------------------------------------------
const JANPACKVec& JANPACKVec::operator*= (const real a)
{
  dolfin_assert(x);
  vec_scal(a, x);
  
  return *this;
}
//-----------------------------------------------------------------------------
const JANPACKVec& JANPACKVec::operator/= (const real a)
{
  dolfin_assert(x);
  dolfin_assert(a != 0.0);
  
  const real b = 1.0 / a;
  vec_scal(b, x);
  
  return *this;
}
//-----------------------------------------------------------------------------
real JANPACKVec::inner(const GenericVector& y) const
{
  dolfin_assert(x);

  const JANPACKVec& v = y.down_cast<JANPACKVec>();
  dolfin_assert(v.x);

  real a;
  a = vec_dot(x, v.x);

  return a;
}
//-----------------------------------------------------------------------------
void JANPACKVec::axpy(real a, const GenericVector& y) 
{
  dolfin_assert(x);

  const JANPACKVec& v = y.down_cast<JANPACKVec>();
  dolfin_assert(v.x);

  vec_axpy(a, v.x, x);
}
//-----------------------------------------------------------------------------
real JANPACKVec::norm(VectorNormType type) const
{
  error("Not implemented.");
  return 0.0;
}
//-----------------------------------------------------------------------------
real JANPACKVec::min() const
{
  error("Not implemented.");
  return 0.0;
}
//-----------------------------------------------------------------------------
real JANPACKVec::max() const
{
  error("Not implemented.");
  return 0.0;
}
//-----------------------------------------------------------------------------
void JANPACKVec::disp(uint precision) const
{
  vec_print(x);
}
//-----------------------------------------------------------------------------  
Vec_ *JANPACKVec::vec() const
{
  return x;
}
//-----------------------------------------------------------------------------  
void JANPACKVec::init_ghosted(uint n, std::set<uint>& indices,
			       std::map<uint, uint>& map)
{
  error("Not implemented.");
}
//-----------------------------------------------------------------------------
LinearAlgebraFactory& JANPACKVec::factory() const
{
  return JANPACKFactory::instance();
}
//-----------------------------------------------------------------------------

#endif
