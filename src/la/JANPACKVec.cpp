// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_JANPACK

#include <dolfin/math/dolfin_math.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/la/JANPACKVec.h>
#include <dolfin/la/JANPACKFactory.h>
#include <dolfin/main/MPI.h>

#include <dolfin/common/Array.h>

#include <cstring>
#include <set>
#include <map>

using namespace dolfin;

//-----------------------------------------------------------------------------
JANPACKVec::JANPACKVec():
    Variable("x", "a sparse vector"),
    is_view(false), is_ghosted(false), is_init(false)
{
  message("Inside vec");
  // Do nothing
}
//-----------------------------------------------------------------------------
JANPACKVec::JANPACKVec(uint N):
    Variable("x", "a sparse vector"), 
    is_view(false), is_ghosted(false), is_init(false)
{
  // Create PETSc vector
  init(N);
}
//-----------------------------------------------------------------------------
JANPACKVec::JANPACKVec(const JANPACKVec& v):
    Variable("x", "a vector"),
    is_view(false), is_ghosted(false), is_init(false)
{
  *this = v;
}
//-----------------------------------------------------------------------------
JANPACKVec::~JANPACKVec()
{
  if (is_init && !is_view)
    jp_vec_free(x);
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
  

  if ((is_init && this->size() == N ) || (this->local_size() == N))
  {
    jp_vec_zero(x);
    return;      
  }
  else
  {
    if (is_init && !is_view) 
    {
      jp_vec_free(x);
    }
    
  }

  // Create vector
  jp_vec_init(x, N);
  is_init = true;
  //  x = &_x;
  jp_vec_zero(x);

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
  //  memcpy(values, x->x, x->n * sizeof(real));
}
//-----------------------------------------------------------------------------
void JANPACKVec::set(real* values)
{
  //  memcpy(x->x,values, x->n * sizeof(real));
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
  jp_vec_get_block(const_cast<char *>(x), const_cast<double*>(block), 
		reinterpret_cast<uint*>(const_cast<uint*>(rows)) , m);
}
//-----------------------------------------------------------------------------
void JANPACKVec::set(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x);
  jp_vec_set_block(const_cast<char *>(x), const_cast<double*>(block), 
		reinterpret_cast<uint*>(const_cast<uint*>(rows)) , m);
}
//-----------------------------------------------------------------------------
void JANPACKVec::add(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x);

  jp_vec_add_block(const_cast<char *>(x), const_cast<double*>(block),
		reinterpret_cast<uint*>(const_cast<uint*>(rows)), m);
  
  //  error("Not implemented.");
}
//-----------------------------------------------------------------------------
void JANPACKVec::apply(FinalizeType finaltype)
{

  jp_vec_finalize(x);
  if (is_ghosted)
    jp_vec_update_ghosts(x);
}
//-----------------------------------------------------------------------------
void JANPACKVec::zero()
{
  dolfin_assert(x);
  jp_vec_zero(x);
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKVec::size() const
{
  uint32_t m = 0;
  jp_vec_size(const_cast<char *>(x), &m);
    return static_cast<uint>(m);
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKVec::local_size() const
{
  uint32_t n = 0;
  jp_vec_local_size(const_cast<char *>(x), &n);

  return static_cast<uint>(n);
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKVec::offset() const
{
  int n = 0;
  //  if (x) 
    //    n = x->range[0];

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
  //jp_vec_copy(x, v.x);
  jp_vec_copy(const_cast<char *>(v.x), const_cast<char *>(x));  
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
  jp_vec_scal(a, x);
  
  return *this;
}
//-----------------------------------------------------------------------------
const JANPACKVec& JANPACKVec::operator/= (const real a)
{
  dolfin_assert(x);
  dolfin_assert(a != 0.0);
  
  const real b = 1.0 / a;
  jp_vec_scal(b, x);
  
  return *this;
}
//-----------------------------------------------------------------------------
real JANPACKVec::inner(const GenericVector& y) const
{
  dolfin_assert(x);

  const JANPACKVec& v = y.down_cast<JANPACKVec>();
  dolfin_assert(v.x);

  real a;
  a = jp_vec_dot(const_cast<char *>(x), const_cast<char *>(v.x));

  return a;
}
//-----------------------------------------------------------------------------
void JANPACKVec::axpy(real a, const GenericVector& y) 
{
  dolfin_assert(x);

  const JANPACKVec& v = y.down_cast<JANPACKVec>();
  dolfin_assert(v.x);

  jp_vec_axpy(a, v.vec(), const_cast<char*>(x));
}
//-----------------------------------------------------------------------------
real JANPACKVec::norm(VectorNormType type) const
{
  return jp_vec_nrm2(const_cast<char*>(x));
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
  jp_vec_print(const_cast<char*>(x));
}
//-----------------------------------------------------------------------------  
char *JANPACKVec::vec() const
{
  return const_cast<char *>(x);
}
//-----------------------------------------------------------------------------  
void JANPACKVec::init_ghosted(uint n, std::set<uint>& indices,
			       std::map<uint, uint>& map)
{
 
  /*
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

  jp_vec_init_ghosts(x, &ghost_indices[0], ghost_indices.size());
  
  is_ghosted = true;
  apply();
  */
}
//-----------------------------------------------------------------------------
LinearAlgebraFactory& JANPACKVec::factory() const
{
  return JANPACKFactory::instance();
}
//-----------------------------------------------------------------------------

#endif
