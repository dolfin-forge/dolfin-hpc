// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_JANPACK

#include <dolfin/math/basic.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/la/JANPACKVec.h>
#include <dolfin/la/JANPACKFactory.h>

#include <dolfin/common/Array.h>

#include <string>

using namespace dolfin;

//-----------------------------------------------------------------------------
JANPACKVec::JANPACKVec():
    Variable("x", "a sparse vector"),
    is_ghosted(false), is_init(false)
#ifdef HAVE_JANPACK_MPI
    , x(&_x)
#endif
{
  // Do nothing
}
//-----------------------------------------------------------------------------
JANPACKVec::JANPACKVec(uint N, bool distributed):
    Variable("x", "a sparse vector"),
    is_ghosted(false), is_init(false)
{
  // Create PETSc vector
  init(N);
}
//-----------------------------------------------------------------------------
JANPACKVec::JANPACKVec(const JANPACKVec& v):
    Variable("x", "a vector"),
    is_ghosted(false), is_init(false)
{
  *this = v;
}
//-----------------------------------------------------------------------------
JANPACKVec::~JANPACKVec()
{
  if (is_init)
    jp_vec_free(x_);
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


  if ((is_init && this->size() == N ) || (is_init && this->local_size() == N))
  {
    jp_vec_zero(x_);
    return;
  }
  else
  {
    if (is_init)
    {
      jp_vec_free(x_);
    }

  }

  // Create vector
  jp_vec_init(x_, N);
  is_init = true;
  //  x = &_x;
  jp_vec_zero(x_);

}
//-----------------------------------------------------------------------------
void JANPACKVec::init(uint N, bool distributed)
{
  init(N);
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
  jp_vec_get_local(const_cast<jp_vec_type *>(x_), values);
}
//-----------------------------------------------------------------------------
void JANPACKVec::set(real* values)
{
  jp_vec_set_local(const_cast<jp_vec_type *>(x_), values);
}
//-----------------------------------------------------------------------------
void JANPACKVec::add(real* values)
{
  dolfin_assert(x_);

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
  dolfin_assert(x_);
  jp_vec_get_block(const_cast<jp_vec_type *>(x_), const_cast<double*>(block),
		   reinterpret_cast<uint*>(const_cast<uint*>(rows)) , m);
}
//-----------------------------------------------------------------------------
void JANPACKVec::set(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x_);
  jp_vec_set_block(const_cast<jp_vec_type *>(x_), const_cast<double*>(block),
		   reinterpret_cast<uint*>(const_cast<uint*>(rows)) , m);
}
//-----------------------------------------------------------------------------
void JANPACKVec::add(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x_);

  jp_vec_add_block(const_cast<jp_vec_type *>(x_), const_cast<double*>(block),
		   reinterpret_cast<uint*>(const_cast<uint*>(rows)), m);

}
//-----------------------------------------------------------------------------
void JANPACKVec::apply(FinalizeType finaltype)
{

  jp_vec_finalize(x_);
  if (is_ghosted)
    jp_vec_update_ghosts(x_);

}
//-----------------------------------------------------------------------------
void JANPACKVec::zero()
{
  dolfin_assert(x_);
  jp_vec_zero(x_);
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKVec::size() const
{
  uint32_t m = 0;
  jp_vec_size(const_cast<jp_vec_type *>(x_), &m);
    return static_cast<uint>(m);
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKVec::local_size() const
{
  uint32_t n = 0;
  jp_vec_local_size(const_cast<jp_vec_type *>(x_), &n);

  return static_cast<uint>(n);
}
//-----------------------------------------------------------------------------
dolfin::uint JANPACKVec::offset() const
{
  uint32_t range[2];

  if(x_)
    jp_vec_range(const_cast<jp_vec_type *>(x_), range);

  return static_cast<uint>(range[0]);
}
//-----------------------------------------------------------------------------
GenericVector& JANPACKVec::operator= (const GenericVector& v)
{
  *this = v.down_cast<JANPACKVec>();
  return *this;
}
//-----------------------------------------------------------------------------
JANPACKVec& JANPACKVec::operator= (const JANPACKVec& v)
{
  dolfin_assert(v.x_);

  init(v.local_size());
  jp_vec_copy(const_cast<jp_vec_type *>(v.x_), const_cast<jp_vec_type *>(x_));
  return *this;
}
//-----------------------------------------------------------------------------
JANPACKVec& JANPACKVec::operator= (real a)
{
  dolfin_assert(x_);
  // VecSet(x, a);
  //  error("Not implemented");
  message("....");
  return *this;
}
//-----------------------------------------------------------------------------
JANPACKVec& JANPACKVec::operator*= (const GenericVector& x)
{
  dolfin_assert(x_);

  const JANPACKVec& v = x.down_cast<JANPACKVec>();
  dolfin_assert(v.x_);

  jp_vec_pwmul(const_cast<jp_vec_type *>(x_),
	       const_cast<jp_vec_type *>(x_), const_cast<jp_vec_type *>(v.x_));

  return *this;
}
//-----------------------------------------------------------------------------
JANPACKVec& JANPACKVec::operator+= (const GenericVector& x)
{
  this->axpy(1.0, x);
  return *this;
}
//-----------------------------------------------------------------------------
JANPACKVec& JANPACKVec::operator-= (const GenericVector& x)
{
  this->axpy(-1.0, x);
  return *this;
}
//-----------------------------------------------------------------------------
JANPACKVec& JANPACKVec::operator*= (const real a)
{
  dolfin_assert(x_);
  jp_vec_scal(a, x_);

  return *this;
}
//-----------------------------------------------------------------------------
JANPACKVec& JANPACKVec::operator/= (const real a)
{
  dolfin_assert(x_);
  dolfin_assert(a != 0.0);

  const real b = 1.0 / a;
  jp_vec_scal(b, x_);

  return *this;
}
//-----------------------------------------------------------------------------
real JANPACKVec::inner(const GenericVector& y) const
{
  dolfin_assert(x_);

  const JANPACKVec& v = y.down_cast<JANPACKVec>();
  dolfin_assert(v.x_);

  real a;
  a = jp_vec_dot(const_cast<jp_vec_type *>(x_),
		 const_cast<jp_vec_type *>(v.x_));

  return a;
}
//-----------------------------------------------------------------------------
void JANPACKVec::axpy(real a, const GenericVector& y)
{
  dolfin_assert(x_);

  const JANPACKVec& v = y.down_cast<JANPACKVec>();
  dolfin_assert(v.x_);

  jp_vec_axpy(a, v.vec(), const_cast<jp_vec_type *>(x_));
}
//-----------------------------------------------------------------------------
real JANPACKVec::norm(VectorNormType type) const
{
  return jp_vec_nrm2(const_cast<jp_vec_type *>(x_));
}
//-----------------------------------------------------------------------------
real JANPACKVec::min() const
{
  return jp_vec_min(const_cast<jp_vec_type *>(x_));
}
//-----------------------------------------------------------------------------
real JANPACKVec::max() const
{
  return jp_vec_max(const_cast<jp_vec_type *>(x_));
}
//-----------------------------------------------------------------------------
void JANPACKVec::pointwise(const GenericVector& x, VectorPointwiseOp op) const
{
  const JANPACKVec& v = x.down_cast<JANPACKVec>();
  dolfin_assert(v.x_);

  switch(op)
  {
  case pw_min:
    jp_vec_pwmin(const_cast<jp_vec_type *>(x_),
		 const_cast<jp_vec_type *>(x_), const_cast<jp_vec_type *>(v.x_));
    break;
  case pw_max:
    jp_vec_pwmax(const_cast<jp_vec_type *>(x_),
		 const_cast<jp_vec_type *>(x_), const_cast<jp_vec_type *>(v.x_));
    break;
  case pw_mult:
    jp_vec_pwmul(const_cast<jp_vec_type *>(x_),
		 const_cast<jp_vec_type *>(x_), const_cast<jp_vec_type *>(v.x_));
    break;
  case pw_div:
    jp_vec_pwdiv(const_cast<jp_vec_type *>(x_),
		 const_cast<jp_vec_type *>(x_), const_cast<jp_vec_type *>(v.x_));
    break;
  default:
    error("Unknown operator");
  }

}
//-----------------------------------------------------------------------------
void JANPACKVec::disp(uint precision) const
{
  jp_vec_print(const_cast<jp_vec_type *>(x_));
}
//-----------------------------------------------------------------------------
jp_vec_type *JANPACKVec::vec() const
{
  return const_cast<jp_vec_type *>(x_);
}
//-----------------------------------------------------------------------------
void JANPACKVec::init_ghosted(uint n, _ordered_set<uint>& indices,
			       _ordered_map<uint, uint>& map)
{

  if ( is_ghosted )
    apply();

  uint32_t range[2];
  jp_vec_range(x_, range);

  Array<uint32_t> ghost_indices;
  _ordered_set<uint32_t>::iterator sit;
  for(sit = indices.begin(); sit != indices.end(); ++sit) {
    if( *sit < (uint32_t) range[0] || *sit >= (uint32_t) range[1] ) {
      ghost_indices.push_back((uint32_t) *sit);
    }
  }

  jp_vec_init_ghosts(x_, &ghost_indices[0], ghost_indices.size());

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
