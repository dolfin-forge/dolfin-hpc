// Copyright (C) 2004-2007 Johan Hoffman, Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2005-2007.
// Modified by Martin Sandve Alnes 2008.
// Modified by Niclas Jansson 2008.
// Modified by Niyazi Cem Degirmenci 2013.
//
// First added:  2004
// Last changed: 2008-10-28

#include <dolfin/la/PETScVector.h>

#ifdef HAVE_PETSC

#include <dolfin/common/Array.h>
#include <dolfin/log/log.h>
#include <dolfin/la/PETScFactory.h>
#include <dolfin/main/PE.h>
#include <dolfin/main/MPI.h>
#include <dolfin/math/basic.h>

#include <set>
#include <map>

namespace dolfin
{

//-----------------------------------------------------------------------------
PETScVector::PETScVector() :
    Variable("x", "a sparse vector"),
    x_(0),
    is_distributed_(false),
    is_ghosted_(false)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(uint N, bool distributed) :
    Variable("x", "a sparse vector"),
    x_(0),
    is_distributed_(false),
    is_ghosted_(false)
{
  // Create PETSc vector
  init(N, distributed);
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(Vec x) :
    Variable("x", "a vector"),
    x_(x),
    is_distributed_(false),
    is_ghosted_(false)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(PETScVector const& v) :
    Variable("x", "a vector"),
    x_(0),
    is_distributed_(false),
    is_ghosted_(false)
{
  *this = v;
}
//-----------------------------------------------------------------------------
PETScVector::~PETScVector()
{
  clear();
}
//-----------------------------------------------------------------------------
void PETScVector::clear()
{
  if (x_)
  {
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    VecDestroy(&x_);
#else
    VecDestroy(x_);
#endif
    is_ghosted_ = false;
    is_distributed_ = false;
  }
}
//-----------------------------------------------------------------------------
void PETScVector::init(uint N)
{
  init(N, true);
}
//-----------------------------------------------------------------------------
void PETScVector::init(uint N, bool distributed)
{
  // Do not reallocate
  if (x_ && this->local_size() == N) { return; }

  clear();

  // Create vector
  if (distributed && PE::size() > 1)
  {
    is_distributed_ = true;
#ifdef HAVE_MPI
    VecCreateMPI(MPI::DOLFIN_COMM, N, PETSC_DETERMINE, &x_);
#endif
  }
  else
  {
    VecCreate(PETSC_COMM_SELF, &x_);
    VecSetSizes(x_, PETSC_DECIDE, N);
    VecSetFromOptions(x_);
  }
  // Set all entries to zero
  PetscScalar a = 0.0;
  VecSet(x_, a);
}
//-----------------------------------------------------------------------------
PETScVector* PETScVector::copy() const
{
  return new PETScVector(*this);
}
//-----------------------------------------------------------------------------
void PETScVector::get(real* values) const
{
#if PETSC_VERSION_MAJOR > 2
  dolfin_assert(x_);
  real const* data = NULL;
  VecGetArrayRead(x_, &data);
  dolfin_assert(data);
  PetscInt n;
  VecGetLocalSize(x_, &n);
  std::copy(data, data + n, values);
  VecRestoreArrayRead(x_, &data);
  dolfin_assert(x_);
#else
  dolfin_assert(x_);

  real* data = 0;
  VecGetArray(x_, &data);
  dolfin_assert(data);

  for (uint i = 0; i < local_size(); i++)
    values[i] = data[i];
  VecRestoreArray(x_, &data);

  dolfin_assert(x_);
#endif
}
//-----------------------------------------------------------------------------
void PETScVector::set(real* values)
{
  dolfin_assert(x_);
  real* data = NULL;
  VecGetArray(x_, &data);
  dolfin_assert(data);
  PetscInt n;
  VecGetLocalSize(x_, &n);
  std::copy(values, values + n, data);
  VecRestoreArray(x_, &data);
  dolfin_assert(x_);
}
//-----------------------------------------------------------------------------
void PETScVector::add(real* values)
{
  dolfin_assert(x_);
  PetscInt n;
  VecGetLocalSize(x_, &n);
  int * rows = new int[n];
  for (int i = 0; i < n; i++) { rows[i] = i; }
  VecSetValues(x_, n, rows, values, ADD_VALUES);
  delete[] rows;
}
//-----------------------------------------------------------------------------
void PETScVector::get(real* block, uint m, const uint* rows) const
{
  dolfin_assert(x_);

  if (is_ghosted_)
  {
    int low, high;
    Vec xl;
    VecGetOwnershipRange(x_, &low, &high);
    VecGhostGetLocalForm(x_, &xl);

    int *tmp = new int[m];
    for (uint i = 0; i < m; i++)
    {
      if ((int) rows[i] < high && (int) rows[i] >= low)
      {
        tmp[i] = rows[i] - low;
      }
      else
      {
        dolfin_assert(mapping_.size() > 0);
        GhostMapping::const_iterator it = mapping_.find(rows[i]);
        dolfin_assert(mapping_.count(rows[i]) > 0);
        tmp[i] = it->second;
      }
    }
    VecGetValues(xl, static_cast<int>(m), tmp, block);
    VecGhostRestoreLocalForm(x_, &xl);

    delete[] tmp;
  }
  else
  {
    VecGetValues(x_, static_cast<int>(m),
                 reinterpret_cast<int*>(const_cast<uint*>(rows)), block);
  }

}
//-----------------------------------------------------------------------------
void PETScVector::set(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x_);
  VecSetValues(x_, static_cast<int>(m),
               reinterpret_cast<int*>(const_cast<uint*>(rows)), block,
               INSERT_VALUES);
}
//-----------------------------------------------------------------------------
void PETScVector::add(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x_);
  VecSetValues(x_, static_cast<int>(m),
               reinterpret_cast<int*>(const_cast<uint*>(rows)), block,
               ADD_VALUES);
}
//-----------------------------------------------------------------------------
void PETScVector::apply(FinalizeType finaltype)
{

  VecAssemblyBegin(x_);
  VecAssemblyEnd(x_);

  if (is_ghosted_)
  {
    VecGhostUpdateBegin(x_, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x_, INSERT_VALUES, SCATTER_FORWARD);
  }
}
//-----------------------------------------------------------------------------
void PETScVector::zero()
{
  dolfin_assert(x_);
  real a = 0.0;
  VecSet(x_, a);
}
//-----------------------------------------------------------------------------
uint PETScVector::size() const
{
  if(x_ == NULL) return 0;
  int n = 0;
  VecGetSize(x_, &n);
  return static_cast<uint>(n);
}
//-----------------------------------------------------------------------------
uint PETScVector::local_size() const
{
  dolfin_assert(x_);
  int n = 0;
  VecGetLocalSize(x_, &n);
  return static_cast<uint>(n);
}
//-----------------------------------------------------------------------------
uint PETScVector::offset() const
{
  dolfin_assert(x_);
  int low, high;
  VecGetOwnershipRange(x_, &low, &high);
  return static_cast<uint>(low);
}
//-----------------------------------------------------------------------------
PETScVector& PETScVector::operator=(const GenericVector& v)
{
  *this = v.down_cast<PETScVector>();
  return *this;
}
//-----------------------------------------------------------------------------
PETScVector& PETScVector::operator=(PETScVector const& v)
{
  if (&v != this)
  {
    dolfin_assert(v.x_);
    init(v.local_size(), v.is_distributed_);
    VecCopy(v.x_, x_);
  }
  return *this;
}
//-----------------------------------------------------------------------------
PETScVector& PETScVector::operator=(real a)
{
  dolfin_assert(x_);
  VecSet(x_, a);
  return *this;
}
//-----------------------------------------------------------------------------
PETScVector& PETScVector::operator*=(const GenericVector& y)
{
  dolfin_assert(x_);
  PETScVector const& v = y.down_cast<PETScVector>();
  dolfin_assert(v.x_);

  if (size() != v.size())
  {
    error("Vectors must have the same size for componentwise multiplication.");
  }

  VecPointwiseMult(x_, x_, v.x_);

  return *this;
}
//-----------------------------------------------------------------------------
PETScVector& PETScVector::operator+=(const GenericVector& x)
{
  this->axpy(1.0, x);
  return *this;
}
//-----------------------------------------------------------------------------
PETScVector& PETScVector::operator-=(const GenericVector& x)
{
  this->axpy(-1.0, x);
  return *this;
}
//-----------------------------------------------------------------------------
PETScVector& PETScVector::operator*=(const real a)
{
  dolfin_assert(x_);
  VecScale(x_, a);

  return *this;
}
//-----------------------------------------------------------------------------
PETScVector& PETScVector::operator/=(const real a)
{
  dolfin_assert(x_);
  dolfin_assert(a != 0.0);

  const real b = 1.0 / a;
  VecScale(x_, b);

  return *this;
}
//-----------------------------------------------------------------------------
real PETScVector::inner(const GenericVector& y) const
{
  dolfin_assert(x_);

  PETScVector const& v = y.down_cast<PETScVector>();
  dolfin_assert(v.x_);

  real a;
  VecDot(v.x_, x_, &a);

  return a;
}
//-----------------------------------------------------------------------------
void PETScVector::axpy(real a, const GenericVector& y)
{
  dolfin_assert(x_);

  PETScVector const& v = y.down_cast<PETScVector>();
  dolfin_assert(v.x_);

  if (size() != v.size())
  {
    error("The vectors must be of the same size to apply AXPY.");
  }

  VecAXPY(x_, a, v.x_);
}
//-----------------------------------------------------------------------------
real PETScVector::norm(VectorNormType type) const
{
  dolfin_assert(x_);

  real value = 0.0;

  switch (type)
    {
    case l1:
      VecNorm(x_, NORM_1, &value);
      break;
    case l2:
      VecNorm(x_, NORM_2, &value);
      break;
    default:
      VecNorm(x_, NORM_INFINITY, &value);
      break;
    }

  return value;
}
//-----------------------------------------------------------------------------
real PETScVector::min() const
{
  real value = 0.0;

  VecMin(x_, PETSC_NULL, &value);

  return value;
}
//-----------------------------------------------------------------------------
real PETScVector::max() const
{
  real value = 0.0;

  VecMax(x_, PETSC_NULL, &value);

  return value;
}
//-----------------------------------------------------------------------------
void PETScVector::disp(uint precision) const
{
  section("PETScVector");
  if (PE::size() > 1 && is_distributed_)
  {
    VecView(x_, PETSC_VIEWER_STDOUT_WORLD);
  }
  else
  {
    VecView(x_, PETSC_VIEWER_STDOUT_SELF);
  }
  end();
}
//-----------------------------------------------------------------------------
Vec PETScVector::vec() const
{
  return x_;
}
//-----------------------------------------------------------------------------
void PETScVector::init_ghosted(uint n, std::set<uint>& indices,
                               std::map<uint, uint>& map)
{
  if (!is_distributed_)
  {
    return;
  }

  if (is_ghosted_)
  {
    apply();
  }

  int local_size, size, low, high;
  VecGetSize(x_, &size);
  VecGetLocalSize(x_, &local_size);
  VecGetOwnershipRange(x_, &low, &high);

  mapping_.clear();

  int *rows = new int[local_size];
  real *values = new real[local_size];
  for (int i = 0; i < local_size; i++)
  {
    rows[i] = low + i;
    mapping_[low + i] = i;
  }

  VecGetValues(x_, local_size, rows, values);

  if (is_ghosted_ && map.size() > 0)
  {
    // dolfin_assert(map.size() > 0);
    for (int i = 0; i < local_size; i++)
    {
      // dolfin_assert(map.count(low + i) > 0);
      rows[i] = map[low + i];
    }
  }

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  VecDestroy(&x_);
#else
  VecDestroy(x_);
#endif

  Array<int> ghost_indices;
  int num_ghost = local_size;
  std::set<uint>::iterator sit;
  for (sit = indices.begin(); sit != indices.end(); ++sit)
  {
    if (*sit < (uint) low || *sit >= (uint) high)
    {
      ghost_indices.push_back((int) *sit);
      mapping_[(int) *sit] = num_ghost++;
    }
  }

#ifdef HAVE_MPI
  VecCreateGhost(MPI::DOLFIN_COMM, local_size, size, (int) ghost_indices.size(),
                 (const int *) &ghost_indices[0], &x_);
#endif
  VecSetValues(x_, local_size, rows, values, INSERT_VALUES);
  delete[] rows;
  delete[] values;

  is_ghosted_ = true;
  apply();
}
//-----------------------------------------------------------------------------
LinearAlgebraFactory& PETScVector::factory() const
{
  return PETScFactory::instance();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* HAVE_PETSC */

