// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-19
// Last changed: 2006-10-19

#include <dolfin/mesh/MeshGeometry.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/log/dolfin_log.h>

#include <cstring>
#include <ctime>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry() :
    dim_(0),
    size_(0),
    coordinates_(NULL),
    timestamp_(std::time(NULL))
{
  std::memset(abs_tol_,0,sizeof(abs_tol_));
}
//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry(uint gdim, uint size) :
    dim_(gdim),
    size_(size),
    coordinates_(NULL),
    timestamp_(std::time(NULL))
{
  std::memset(abs_tol_, 0, sizeof(abs_tol_));
  init(gdim, size);
}
//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry(MeshGeometry const& geometry) :
    dim_(0),
    size_(0),
    coordinates_(NULL),
    timestamp_(std::time(NULL))
{
  *this = geometry;
}
//-----------------------------------------------------------------------------
MeshGeometry::~MeshGeometry()
{
  clear();
}
//-----------------------------------------------------------------------------
MeshGeometry const& MeshGeometry::operator=(MeshGeometry const& geometry)
{
  // Clear old data if any
  clear();

  // Allocate data
  dim_ = geometry.dim_;
  size_ = geometry.size_;
  uint const n = dim_ * size_;
  coordinates_ = new real[n];

  // Copy data
  for (uint i = 0; i < n; ++i)
  {
    coordinates_[i] = geometry.coordinates_[i];
  }

  // Copy tolerances
  for (uint i = 0; i <= dim_; ++i)
  {
    abs_tol_[i] = geometry.abs_tol_[i];
  }

  timestamp_ = geometry.timestamp_;

  return *this;
}
//-----------------------------------------------------------------------------
uint MeshGeometry::dim() const
{
  return dim_;
}
//-----------------------------------------------------------------------------
uint MeshGeometry::size() const
{
  return size_;
}
//-----------------------------------------------------------------------------
real MeshGeometry::abs_tolerance(uint dim) const
{
  dolfin_assert(dim <= dim_);
  return abs_tol_[dim];
}
//-----------------------------------------------------------------------------
Point MeshGeometry::point(uint n) const
{
  Point p;
  std::memcpy(&p[0], &coordinates_[n * dim_], dim_*sizeof(real));
  return p;
}
//-----------------------------------------------------------------------------
real * MeshGeometry::coordinates()
{
  return coordinates_;
}
//-----------------------------------------------------------------------------
real const * MeshGeometry::coordinates() const
{
  return coordinates_;
}
//-----------------------------------------------------------------------------
void MeshGeometry::clear()
{
  dim_ = 0;
  size_ = 0;
  delete[] coordinates_;
  std::memset(abs_tol_,0,sizeof(abs_tol_));
  coordinates_ = NULL;
}
//-----------------------------------------------------------------------------
void MeshGeometry::finalize()
{
  dolfin_assert(coordinates_ != NULL);
}
//-----------------------------------------------------------------------------
void MeshGeometry::init(uint gdim, uint size)
{
  // Delete old data if any
  clear();

  // Allocate new data
  coordinates_ = new real[gdim * size];
  std::memset(coordinates_, 0.0, gdim * size * sizeof(real));

  // Save dimension and size
  dim_ = gdim;
  size_ = size;

  // Initialize tolerances
  for(uint i = 0; i <= dim_; ++i)
  {
    abs_tol_[i] = DOLFIN_EPS;
  }

  // Initialize token
  update_token();
}
//-----------------------------------------------------------------------------
void MeshGeometry::set_abs_tolerance(uint dim, real atol)
{
  dolfin_assert(dim <= dim_);
  abs_tol_[dim] = std::fabs(atol);
}
//-----------------------------------------------------------------------------
void MeshGeometry::set(uint n, uint i, real x)
{
  coordinates_[n * dim_ + i] = x;
}
//-----------------------------------------------------------------------------
void MeshGeometry::set(uint n, real const * x)
{
  std::memcpy(&coordinates_[n * dim_], x, dim_*sizeof(real));
}
//-----------------------------------------------------------------------------
void MeshGeometry::remap(Array<uint> const& map)
{
  real * xcpy = new real[size_ * dim_];
  for (uint i = 0; i < size_; ++i)
  {
    std::memcpy(&xcpy[map[i] * dim_], &coordinates_[i * dim_],
                dim_ * sizeof(real));
  }
  delete[] coordinates_;
  coordinates_ = xcpy;

  // Invalidate dependencies
  update_token();
}
//-----------------------------------------------------------------------------
int MeshGeometry::token() const
{
  return timestamp_ + size(); // FIXME
}
//-----------------------------------------------------------------------------
void MeshGeometry::update_token()
{
  timestamp_ = std::time(NULL);
}
//-----------------------------------------------------------------------------
void MeshGeometry::disp() const
{
  // Begin indentation
  cout << "Mesh geometry" << endl;
  begin("-------------");
  cout << endl;

  // Check if empty
  if (dim_ == 0)
  {
    cout << "empty" << endl << endl;
    end();
    return;
  }

  // Display coordinates for all vertices
  for (uint i = 0; i < size_; ++i)
  {
    cout << i << ":";
    for (uint d = 0; d < dim_; ++d)
    {
      cout << " " << x(i, d);
    }
    cout << endl;
  }
  cout << endl;

  // End indentation
  end();
}
//-----------------------------------------------------------------------------
void MeshGeometry::check() const
{
  //FIXME
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

