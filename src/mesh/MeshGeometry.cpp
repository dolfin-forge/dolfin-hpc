// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-19
// Last changed: 2006-10-19

#include <dolfin/mesh/MeshGeometry.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/log/log.h>

#include <cstring>
#include <ctime>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry() :
    dim_(0),
    size_(0),
    coordinates_(NULL),
    abs_tol_(NULL),
    timestamp_(0)
{
}
//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry(uint gdim, uint size) :
    dim_(0),
    size_(0),
    coordinates_(NULL),
    abs_tol_(NULL),
    timestamp_(0)
{
  init(gdim, size);
}
//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry(MeshGeometry const& geometry) :
    dim_(0),
    size_(0),
    coordinates_(NULL),
    timestamp_(0)
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
  clear();

  dim_ = geometry.dim_;
  size_ = geometry.size_;
  if(dim_*size_ > 0)
  {
    coordinates_ = new real[dim_*size_];
    std::memcpy(coordinates_, geometry.coordinates_, dim_*size_*sizeof(real));
  }
  abs_tol_  = new real[dim_+1];
  std::memcpy(abs_tol_, geometry.abs_tol_, (dim_+1)*sizeof(real));
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
  return Point(&coordinates_[n * dim_], dim_);
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
void MeshGeometry::init(uint gdim, uint size)
{
  if (coordinates_ != NULL)
  {
    error("MeshGeometry : clear instance before reinitializing");
  }
  if (gdim == 0)
  {
    error("MeshGeometry : geometric dimension is zero");
  }
  if (gdim > Point::MAX_SIZE)
  {
    error("MeshGeometry : geometric dimension '%u' exceeds point size", gdim);
  }

  dim_ = gdim;
  size_ = size;
  if(dim_*size_ > 0)
  {
    coordinates_ = new real[dim_*size_];
    std::memset(coordinates_, 0.0, dim_* size_*sizeof(real));
  }
  abs_tol_  = new real[dim_+1];
  std::memset(abs_tol_, DOLFIN_EPS, (dim_+1)*sizeof(real));

  // Initialize token
  update_token();
}
//-----------------------------------------------------------------------------
void MeshGeometry::clear()
{
  dim_ = 0;
  size_ = 0;
  delete[] coordinates_;
  coordinates_ = NULL;
  delete[] abs_tol_;
  abs_tol_ = NULL;
}
//-----------------------------------------------------------------------------
void MeshGeometry::finalize()
{
  if(size_ > 0 && (coordinates_ == NULL))
  {
    error("MeshGeometry : empty coordinates for non-empty geometry");
  }
  // Invalidate dependencies
  update_token();
}
//-----------------------------------------------------------------------------
void MeshGeometry::set_abs_tolerance(uint dim, real atol)
{
  dolfin_assert(dim <= dim_);
  if(atol <= 0.0)
  {
    warning("MeshGeometry : tolerance '%g' for dimension %u is non-positive",
            atol, dim);
  }
  abs_tol_[dim] = std::fabs(atol);
}
//-----------------------------------------------------------------------------
void MeshGeometry::set(uint n, uint i, real x)
{
  dolfin_assert(n < size_);
  dolfin_assert(i < dim_);
  coordinates_[n * dim_ + i] = x;
}
//-----------------------------------------------------------------------------
void MeshGeometry::set(uint n, real const * x)
{
  dolfin_assert(n < size_);
  std::memcpy(&coordinates_[n * dim_], x, dim_*sizeof(real));
}
//-----------------------------------------------------------------------------
void MeshGeometry::remap(Array<uint> const& map)
{
  if (map.size() != size_)
  {
    error("MeshGeometry : size mismatch for remapping of coordinates ");
  }

  real * xcpy = new real[dim_*size_];
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
  return timestamp_ ^ size_;
}
//-----------------------------------------------------------------------------
void MeshGeometry::update_token()
{
  timestamp_ = std::time(NULL);
}
//-----------------------------------------------------------------------------
void MeshGeometry::disp() const
{
  section("MeshGeometry");
  //---
  cout << "dimension   : " << dim_ << endl;
  cout << "coordinates : " << endl << endl;
  if (size_ == 0)
  {
    cout << "empty" << endl << endl;
  }
  else
  {
    for (uint i = 0; i < size_; ++i)
    {
      cout << i << ":";
      for (uint d = 0; d < dim_; ++d)
      {
        cout << " " << x(i, d);
      }
      cout << endl;
    }
  }
  //---
  end();
  skip();
}
//-----------------------------------------------------------------------------
void MeshGeometry::check() const
{
  //FIXME: To be implemented
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

