// Copyright (C) 2004-2007 Johan Hoffman, Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2005-2007.
// Modified by Martin Sandve Alnes 2008
// Modified by Niclas Jansson 2008
//
// First added:  2004
// Last changed: 2008-07-03

// FIXME: Insert dolfin_assert() where appropriate

#ifdef HAS_PETSC

#include <cmath>
#include <dolfin/math/dolfin_math.h>
#include <dolfin/log/dolfin_log.h>
#include "PETScVector.h"
#include "uBlasVector.h"
#include "PETScFactory.h"
#include <dolfin/main/MPI.h>

#include <dolfin/common/Array.h>

#include <set>
#include <map>

using namespace dolfin;

//-----------------------------------------------------------------------------
PETScVector::PETScVector():
    Variable("x", "a sparse vector"),
    x(0), is_view(false), is_ghosted(false)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(uint N):
    Variable("x", "a sparse vector"), 
    x(0), is_view(false), is_ghosted(false)
{
  // Create PETSc vector
  init(N);
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(Vec x):
    Variable("x", "a vector"),
    x(x), is_view(true), is_ghosted(false)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(const PETScVector& v):
    Variable("x", "a vector"),
    x(0), is_view(false), is_ghosted(false)
{
  *this = v;
}
//-----------------------------------------------------------------------------
PETScVector::~PETScVector()
{
  if (x && !is_view)
    VecDestroy(x);
}
//-----------------------------------------------------------------------------
void PETScVector::init(uint N)
{
  // Two cases:
  //
  //   1. Already allocated and dimension changes -> reallocate
  //   2. Not allocated -> allocate
  //
  // Otherwise do nothing
  
  if ((x && this->size() == N ) || (this->local_size() == N))
  {
    //    VecZeroEntries(x);
    return;      
  }
  else
  {
    if (x && !is_view)
      VecDestroy(x);
  }

  // Create vector
  if (MPI::numProcesses() > 1)
  {
    //    VecCreateMPI(MPI::DOLFIN_COMM, PETSC_DECIDE, N, &x);
    VecCreateMPI(MPI::DOLFIN_COMM, N, PETSC_DETERMINE, &x);
  }
  else {
    VecCreate(PETSC_COMM_SELF, &x);

    VecSetSizes(x, PETSC_DECIDE, N);
    VecSetFromOptions(x);

  }
  // Set all entries to zero
  PetscScalar a = 0.0;
  VecSet(x, a);
}
//-----------------------------------------------------------------------------
PETScVector* PETScVector::copy() const
{
  PETScVector* v = new PETScVector(*this); 
  return v; 
}
//-----------------------------------------------------------------------------
void PETScVector::get(real* values) const
{

  real* data = 0;
  VecGetArray(x, &data);
  dolfin_assert(data);

  for (uint i = 0; i < local_size(); i++)
    values[i] = data[i];
  VecRestoreArray(x, &data); 

  dolfin_assert(x);

  /*
  int m = static_cast<int>(local_size());
  int* rows = new int[m];
  for (int i = 0; i < m; i++)
    rows[i] = i;

  VecGetValues(x, m, rows, values);

  delete [] rows;
  */
}
//-----------------------------------------------------------------------------
void PETScVector::set(real* values)
{
  real* data = 0;
  VecGetArray(x, &data);
  dolfin_assert(data);
  
  for (uint i = 0; i < local_size(); i++)
    data[i] = values[i];
  VecRestoreArray(x, &data); 

  dolfin_assert(x);

  /*
  dolfin_assert(x);
  
  int m = static_cast<int>(size());
  int* rows = new int[m];
  for (int i = 0; i < m; i++)
    rows[i] = i;

  VecSetValues(x, m, rows, values, INSERT_VALUES);

  delete [] rows;
  */
}
//-----------------------------------------------------------------------------
void PETScVector::add(real* values)
{
  dolfin_assert(x);
  
  int m = static_cast<int>(size());
  int* rows = new int[m];
  for (int i = 0; i < m; i++)
    rows[i] = i;

  VecSetValues(x, m, rows, values, ADD_VALUES);

  delete [] rows;
}
//-----------------------------------------------------------------------------
void PETScVector::get(real* block, uint m, const uint* rows) const
{
  dolfin_assert(x);

  if( is_ghosted ) {    
    int  low, high;
    Vec xl;
    VecGetOwnershipRange(x, &low, &high);
    VecGhostGetLocalForm(x, &xl);

    int *tmp = new int[m];
    for(uint i = 0; i < m; i++)
      if( (int) rows[i] < high && (int) rows[i] >= low) 
	tmp[i] = rows[i] - low;
      else  {
	dolfin_assert(mapping.size() > 0);
#if (sun || __sun)
	std::map<int, int>::const_iterator it = mapping.find(rows[i]);    
#else
	std::map<const int, int>::const_iterator it = mapping.find(rows[i]);    
#endif
	dolfin_assert(mapping.count(rows[i]) > 0);
	tmp[i] = it->second;
      }
    
    VecGetValues(xl, static_cast<int>(m), tmp, block);
    VecGhostRestoreLocalForm(x, &xl);
    
    delete[] tmp;   
  }
  else
    VecGetValues(x, static_cast<int>(m), reinterpret_cast<int*>(const_cast<uint*>(rows)), block);

}
//-----------------------------------------------------------------------------
void PETScVector::set(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x);
  VecSetValues(x, static_cast<int>(m), reinterpret_cast<int*>(const_cast<uint*>(rows)), block,
	       INSERT_VALUES);
}
//-----------------------------------------------------------------------------
void PETScVector::add(const real* block, uint m, const uint* rows)
{
  dolfin_assert(x);
  VecSetValues(x, static_cast<int>(m), reinterpret_cast<int*>(const_cast<uint*>(rows)), block,
               ADD_VALUES);
}
//-----------------------------------------------------------------------------
void PETScVector::apply(FinalizeType finaltype)
{

  VecAssemblyBegin(x);
  VecAssemblyEnd(x);

  if( is_ghosted ) {
    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
  }
}
//-----------------------------------------------------------------------------
void PETScVector::zero()
{
  dolfin_assert(x);
  real a = 0.0;
  VecSet(x, a);
}
//-----------------------------------------------------------------------------
dolfin::uint PETScVector::size() const
{
  int n = 0;
  if (x)
    VecGetSize(x, &n);
  
  return static_cast<uint>(n);
}
//-----------------------------------------------------------------------------
dolfin::uint PETScVector::local_size() const
{
  int n = 0;
  if (x) 
    VecGetLocalSize(x, &n);
  
  return static_cast<uint>(n);
}
//-----------------------------------------------------------------------------
const GenericVector& PETScVector::operator= (const GenericVector& v)
{
  *this = v.down_cast<PETScVector>();
  return *this; 
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator= (const PETScVector& v)
{
  dolfin_assert(v.x);

  init(v.local_size());
  VecCopy(v.x, x);

  return *this; 
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator= (real a)
{
  dolfin_assert(x);
  VecSet(x, a);
  return *this; 
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator+= (const GenericVector& x)
{
  this->axpy(1.0, x); 
  return *this;
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator-= (const GenericVector& x)
{
  this->axpy(-1.0, x); 
  return *this;
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator*= (const real a)
{
  dolfin_assert(x);
  VecScale(x, a);
  
  return *this;
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator/= (const real a)
{
  dolfin_assert(x);
  dolfin_assert(a != 0.0);

  const real b = 1.0 / a;
  VecScale(x, b);
  
  return *this;
}
//-----------------------------------------------------------------------------
real PETScVector::inner(const GenericVector& y) const
{
  dolfin_assert(x);

  const PETScVector& v = y.down_cast<PETScVector>();
  dolfin_assert(v.x);

  real a;
  VecDot(v.x, x, &a);

  return a;
}
//-----------------------------------------------------------------------------
void PETScVector::axpy(real a, const GenericVector& y) 
{
  dolfin_assert(x);

  const PETScVector& v = y.down_cast<PETScVector>();
  dolfin_assert(v.x);

  if (size() != v.size())
    error("The vectors must be of the same size.");  

  VecAXPY(x, a, v.x);
}
//-----------------------------------------------------------------------------
real PETScVector::norm(VectorNormType type) const
{
  dolfin_assert(x);

  real value = 0.0;

  switch (type) {
  case l1:
    VecNorm(x, NORM_1, &value);
    break;
  case l2:
    VecNorm(x, NORM_2, &value);
    break;
  default:
    VecNorm(x, NORM_INFINITY, &value);
  }
  
  return value;
}
//-----------------------------------------------------------------------------
real PETScVector::min() const
{
  real value = 0.0;

  VecMin(x, PETSC_NULL, &value);

  return value;
}
//-----------------------------------------------------------------------------
real PETScVector::max() const
{
  real value = 0.0;

  VecMax(x, PETSC_NULL, &value);

  return value;
}
//-----------------------------------------------------------------------------
void PETScVector::disp(uint precision) const
{
  if(MPI::numProcesses() > 1)
    VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  else
    VecView(x, PETSC_VIEWER_STDOUT_SELF);
}
//-----------------------------------------------------------------------------  
Vec PETScVector::vec() const
{
  return x;
}
//-----------------------------------------------------------------------------
void PETScVector::init_ghosted(uint n, std::set<uint>& indices,
			       std::map<uint, uint>& map){
 
  if( is_ghosted )
    apply();

  int local_size, size, low, high;
  VecGetSize(x, &size);
  VecGetLocalSize(x, &local_size);
  VecGetOwnershipRange(x, &low, &high);

  mapping.clear();

  int *rows = new int[ local_size ];
  real *values = new real[ local_size ]; 
  for(int i = 0; i < local_size; i++)   {
    rows[i] = low + i;
    mapping[ low + i ] = i;
  }

  VecGetValues(x, local_size, rows, values);    
  
  if( is_ghosted && map.size() > 0) {
    //    dolfin_assert(map.size() > 0);
    for(int i = 0; i < local_size; i++)  {
      // dolfin_assert(map.count(low + i) > 0);
      rows[i] = map[low + i];
    }
  }

  VecDestroy(x);
  
  Array<int> ghost_indices;
  int num_ghost = local_size;
  std::set<uint>::iterator sit;
  for(sit = indices.begin(); sit != indices.end(); ++sit) {
    if( *sit < (uint) low || *sit >= (uint) high ) {
      ghost_indices.push_back((int) *sit);
      mapping[ (int) *sit ] = num_ghost++;
    }
  }

  VecCreateGhost(MPI::DOLFIN_COMM, local_size, size, (int) ghost_indices.size(),
		 (const int *) &ghost_indices[0], &x);       
  VecSetValues(x, local_size, rows, values, INSERT_VALUES);
  delete[] rows;
  delete[] values;

  is_ghosted = true;
  apply();
}
//-----------------------------------------------------------------------------
LinearAlgebraFactory& PETScVector::factory() const
{
  return PETScFactory::instance();
}
//-----------------------------------------------------------------------------
LogStream& dolfin::operator<< (LogStream& stream, const PETScVector& x)
{
  stream << "[ PETSc vector of size " << x.size() << " ]";
  return stream;
}
//-----------------------------------------------------------------------------


#endif
