// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-19
// Last changed: 2006-10-19

#include <dolfin/mesh/MeshGeometry.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/Array.h>
#include <dolfin/math/basic.h>
#include <dolfin/log/log.h>

#include <algorithm>
#include <cstring>
#include <ctime>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry() :
    space_(NULL),
    dim_(0),
    size_(0),
    coordinates_(NULL),
    abs_tol_(NULL),
    timestamp_(0)
{
}
//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry(MeshGeometry const& geometry) :
    space_(NULL),
    dim_(0),
    size_(0),
    coordinates_(NULL),
    abs_tol_(NULL),
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
MeshGeometry const& MeshGeometry::operator=(MeshGeometry const& other)
{
  clear();

  space_ = other.space_->clone();
  dim_ = other.dim_;
  size_ = other.size_;
  uint const n = dim_ * size_;
  if (n > 0)
  {
    coordinates_ = new real[n];
    std::copy(other.coordinates_, other.coordinates_ + n, coordinates_);
  }
  if (other.abs_tol_ != NULL)
  {
    abs_tol_  = new real[dim_+1];
    std::copy(other.abs_tol_, other.abs_tol_ + (dim_ + 1), abs_tol_);
  }
  timestamp_ = other.timestamp_;

  return *this;
}
//-----------------------------------------------------------------------------
bool MeshGeometry::operator==(MeshGeometry const& other) const
{
  if (!objptrcmp(space_, other.space_))
  {
    return false;
  }
  if (size_ != other.size_)
  {
    return false;
  }
  return cmp<real>(dim_ * size_, coordinates_, other.coordinates_);
}
//-----------------------------------------------------------------------------
bool MeshGeometry::operator!=(MeshGeometry const& other) const
{
  return !(*this == other);
}
//-----------------------------------------------------------------------------
void MeshGeometry::swap(MeshGeometry& other)
{
  if (this == &other) return;
  std::swap(space_, other.space_);
  std::swap(dim_, other.dim_);
  std::swap(size_, other.size_);
  std::swap(coordinates_, other.coordinates_);
  std::swap(abs_tol_, other.abs_tol_);
  std::swap(timestamp_, other.timestamp_);
}
//-----------------------------------------------------------------------------
Space const& MeshGeometry::space() const
{
  dolfin_assert(space_);
  return *space_;
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
  return Point(dim_, coordinates_+ n * dim_);
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
void MeshGeometry::init(Space const& space, uint size)
{
  if (coordinates_ != NULL)
  {
    error("MeshGeometry : clear instance before reinitializing");
  }
  space_ = space.clone();
  dim_ = space.dim();
  if (dim_ == 0)
  {
    error("MeshGeometry : geometric dimension is zero");
  }
  if (dim_ > Point::MAX_SIZE)
  {
    error("MeshGeometry : geometric dimension '%u' exceeds point size", dim_);
  }
  size_ = size;
  if (dim_ * size_ > 0)
  {
    coordinates_ = new real[dim_ * size_];
    std::fill_n(coordinates_, dim_ * size_, 0.0);
  }
  abs_tol_ = new real[dim_ + 1];
  std::fill_n(abs_tol_, (dim_ + 1), DOLFIN_EPS);

  // Initialize token
  update_token();
}
//-----------------------------------------------------------------------------
void MeshGeometry::clear()
{
  delete space_;
  space_ = NULL;
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
void MeshGeometry::assign(real const * x)
{
  std::copy(x, x + dim_ * size_, coordinates_);
}

//-----------------------------------------------------------------------------
void MeshGeometry::assign(Array<real> const& coordinates)
{
  if(coordinates.size() != dim_ * size_)
  {
    error("MeshGeometry : size mismatch in coordinates assignment");
  }
  std::copy(coordinates.begin(), coordinates.end(), coordinates_);
}
//-----------------------------------------------------------------------------
void MeshGeometry::remap(Array<uint> const& mapping)
{
  if (mapping.size() != size_)
  {
    error("MeshGeometry : size mismatch for remapping of coordinates ");
  }

  real * xcpy = new real[dim_*size_];
  for (uint i = 0; i < size_; ++i)
  {
    std::copy(coordinates_ + i * dim_, coordinates_ + (i + 1) * dim_,
              xcpy + mapping[i] * dim_);
  }
  delete[] coordinates_;
  coordinates_ = xcpy;

  // Invalidate dependencies
  update_token();
}
//-----------------------------------------------------------------------------
void MeshGeometry::assign(MeshGeometry const& other, Array<uint> const& mapping)
{
  if (this == &other)
  {
    error("MeshGeometry : assignment of coordinates to self");
  }
  if (other.dim() != dim_)
  {
    error("MeshGeometry : dimension mismatch for assignment of coordinates ");
  }
  if (mapping.size() != size_)
  {
    error("MeshGeometry : size mismatch for assignment of coordinates ");
  }
  for (uint i = 0; i < mapping.size(); ++i)
  {
    std::copy(other.x(mapping[i]), other.x(mapping[i]) + dim_,
              coordinates_ + i * dim_);
  }
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator*=(real const a)
{
  real * it  = coordinates_;
  real const * const end = coordinates_ + dim_ * size_;
  while (it != end)
  {
    (*it++) *= a;
  }
  return *this;
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator/=(real const a)
{
  real * it  = coordinates_;
  real const * const end = coordinates_ + dim_ * size_;
  if(small(a))
  {
    error("MeshGeometry : dividing coordinates by zero");
  }
  real const b = 1.0/a;
  while (it != end)
  {
    (*it++) *= b;
  }
  return *this;
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator+=(real const a)
{
  real * it  = coordinates_;
  real const * const end = coordinates_ + dim_ * size_;
  while (it != end)
  {
    (*it++) += a;
  }
  return *this;
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator-=(real const a)
{
  real * it  = coordinates_;
  real const * const end = coordinates_ + dim_ * size_;
  while (it != end)
  {
    (*it++) -= a;
  }
  return *this;
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator+=(Point const& p)
{
  real const N = dim_ * size_;
  for (uint i = 0; i < N; ++i)
  {
    coordinates_[i] += p[i % dim_];
  }
  return *this;
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator-=(Point const& p)
{
  real const N = dim_ * size_;
  for (uint i = 0; i < N; ++i)
  {
    coordinates_[i] -= p[i % dim_];
  }
  return *this;
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
        cout << " " << x(i)[d];
      }
      cout << endl;
    }
  }
  //---
  endblock();
}
//-----------------------------------------------------------------------------
MeshGeometry const& MeshGeometry::operator>>(Array<real>& A) const
{
  A.assign(coordinates_, coordinates_ + dim_ * size_);
  A %= dim_;
  return *this;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

