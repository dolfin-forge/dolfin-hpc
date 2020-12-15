// Copyright (C) 2004-2007 Johan Hoffman, Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/la/PETScVector.h>

#ifdef HAVE_PETSC

#include <dolfin/common/Array.h>
#include <dolfin/la/PETScFactory.h>
#include <dolfin/log/log.h>
#include <dolfin/main/PE.h>
#include <dolfin/main/MPI.h>
#include <dolfin/math/basic.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
PETScVector::PETScVector() :
    Variable("x", "a sparse vector")
    
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(uint N, bool distributed) :
    Variable("x", "a sparse vector"),
    x_(nullptr),
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
    x_(nullptr),
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
void PETScVector::init(uint N, bool distributed)
{
  // Do not reallocate
  if (x_ && this->local_size() == N) { return; }

  clear();

  // Create vector
  if (distributed && PE::size() > 1)
  {
    is_distributed_ = true;
#ifdef DOLFIN_HAVE_MPI
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
void PETScVector::get(real* values) const
{
#if PETSC_VERSION_MAJOR > 2
  dolfin_assert(x_);
  real const* data = nullptr;
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
  real* data = nullptr;
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
void PETScVector::pointwise(const GenericVector& x, VectorPointwiseOp op) const
{

  const PETScVector& v = x.down_cast<PETScVector>();
  dolfin_assert(v.x_);

  switch(op)
  {
  case pw_min:
    VecPointwiseMin(x_, x_, v.x_);
    break;
  case pw_max:
    VecPointwiseMax(x_, x_, v.x_);
    break;
  case pw_mult:
    VecPointwiseMult(x_, x_, v.x_);
    break;
  case pw_div:
    VecPointwiseDivide(x_, x_, v.x_);
    break;
  case pw_maxabs:
    VecPointwiseMaxAbs(x_, x_, v.x_);
    break;
  default:
    error("Unknown operator");
  }

}
//-----------------------------------------------------------------------------
void PETScVector::disp(uint) const
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
void PETScVector::init_ghosted(uint, _ordered_set<uint>& indices,
                               _ordered_map<uint, uint>& map)
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
  _ordered_set<uint>::iterator sit;
  for (sit = indices.begin(); sit != indices.end(); ++sit)
  {
    if (*sit < (uint) low || *sit >= (uint) high)
    {
      ghost_indices.push_back((int) *sit);
      mapping_[(int) *sit] = num_ghost++;
    }
  }

#ifdef DOLFIN_HAVE_MPI
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

