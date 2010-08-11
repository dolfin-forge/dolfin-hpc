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
  zero();
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
  vec_get_block(x, const_cast<double*>(block), 
		reinterpret_cast<int*>(const_cast<uint*>(rows)) , m);
}
//-----------------------------------------------------------------------------
void JANPACKVec::set(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x);
  vec_set_block(x, const_cast<double*>(block), 
		reinterpret_cast<int*>(const_cast<uint*>(rows)) , m);
}
//-----------------------------------------------------------------------------
void JANPACKVec::add(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x);

  vec_add_block(x, const_cast<double*>(block),
		reinterpret_cast<int*>(const_cast<uint*>(rows)), m);
  
  //  error("Not implemented.");
}
//-----------------------------------------------------------------------------
void JANPACKVec::apply(FinalizeType finaltype)
{

  finalize_vec(x);
  if (is_ghosted)
    vec_update_ghosts(x);
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
  //vec_copy(x, v.x);
  vec_copy(v.x, x);
  return *this; 
}
//-----------------------------------------------------------------------------
const JANPACKVec& JANPACKVec::operator= (real a)
{
  dolfin_assert(x);
  // VecSet(x, a);
  error("Not implemented");
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
  return vec_nrm2(x);
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
 
  if ( is_ghosted )
    apply();
  
  int low, high;

  low = x->range[0];
  high = x->range[1];

  Array<int> ghost_indices;
  std::set<uint>::iterator sit;
  for(sit = indices.begin(); sit != indices.end(); ++sit) {
    if( *sit < (uint) low || *sit >= (uint) high ) {
      ghost_indices.push_back((int) *sit);
    }
  }

  init_vec_ghosts(x, &ghost_indices[0], ghost_indices.size());
  
  is_ghosted = true;
  apply();
}
//-----------------------------------------------------------------------------
LinearAlgebraFactory& JANPACKVec::factory() const
{
  return JANPACKFactory::instance();
}
//-----------------------------------------------------------------------------

#endif
