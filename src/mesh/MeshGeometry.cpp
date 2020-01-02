// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

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
MeshGeometry::MeshGeometry(Space const& space, uint size) :
    space_(space.clone()),
    dim_(space.dim()),
    size_(0),
    coordinates_(Array<real>(size_, 0.0)),
    abs_tol_(Array<real>(dim_ + 1, 0.0)),
    timestamp_(0)
{
  resize(size);
}
//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry(MeshGeometry const& other) :
    space_(other.space_->clone()),
    dim_(other.dim_),
    size_(other.size_),
    coordinates_(other.coordinates_),
    abs_tol_(other.abs_tol_),
    timestamp_(other.timestamp_)
{
}
//-----------------------------------------------------------------------------
MeshGeometry::~MeshGeometry()
{
}
//-----------------------------------------------------------------------------
MeshGeometry & MeshGeometry::operator=(MeshGeometry const& other)
{
  MeshGeometry tmp(other);
  swap( *this, tmp );

  return *this;
}
//-----------------------------------------------------------------------------
bool MeshGeometry::operator==(MeshGeometry const& other) const
{
  if ( *space_ != *other.space_ )
  {
    return false;
  }
  if ( size_ != other.size_ )
  {
    return false;
  }
  if ( dim_ != other.dim_ )
  {
    return false;
  }
  if ( coordinates_ != other.coordinates_ )
  {
    return false;
  }
  if ( abs_tol_ != other.abs_tol_ )
  {
    return false;
  }

  return true;
}
//-----------------------------------------------------------------------------
bool MeshGeometry::operator!=(MeshGeometry const& other) const
{
  return !(*this == other);
}
//-----------------------------------------------------------------------------
void swap( MeshGeometry& a, MeshGeometry& b )
{
  using std::swap;

  swap(a.space_       , b.space_);
  swap(a.dim_         , b.dim_);
  swap(a.size_        , b.size_);
  swap(a.coordinates_ , b.coordinates_);
  swap(a.abs_tol_     , b.abs_tol_);
  swap(a.timestamp_   , b.timestamp_);
}
//-----------------------------------------------------------------------------
Space const& MeshGeometry::space() const
{
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
  return Point(dim_, coordinates_.data() + n * dim_);
}
//-----------------------------------------------------------------------------
real * MeshGeometry::coordinates()
{
  return coordinates_.data();
}
//-----------------------------------------------------------------------------
real const * MeshGeometry::coordinates() const
{
  return coordinates_.data();
}
//-----------------------------------------------------------------------------
void MeshGeometry::resize(uint size)
{
  if (size != size_)
  {
    if (size)
    {
      // PointCells have dim == 0, so special treatment is necessary
      coordinates_.resize( ( dim_ > 0 ) ? dim_ * size : size );
      size_ = size;
    }
    else
    {
      coordinates_.clear();
      size_ = 0;
    }
  }
  update_token();
}
//-----------------------------------------------------------------------------
void MeshGeometry::finalize()
{
  if( size_ > 0 && coordinates_.empty() )
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
void MeshGeometry::assign(Array<real> const& coordinates)
{
  if(coordinates.size() % dim_)
  {
    error("MeshGeometry : size mismatch in coordinates assignment");
  }
  coordinates_.resize( coordinates.size() );
  std::copy(coordinates.begin(), coordinates.end(), coordinates_.begin());
  size_ = coordinates.size() / dim_;
}
//-----------------------------------------------------------------------------
void MeshGeometry::remap(Array<uint> const& mapping)
{
  if (mapping.size() != size_)
  {
    error("MeshGeometry : size mismatch for remapping of coordinates ");
  }

  // Reorder coordinates w.r.t old -> new index mapping
  Array<real> xcpy( dim_ * size_, 0.0 );
  for (uint i = 0; i < size_; ++i)
  {
    real const * x = coordinates_.data() + i * dim_;
    std::copy(x, x + dim_, xcpy.data() + mapping[i] * dim_);
  }
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
              coordinates_.data() + i * dim_);
  }
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator*=(real const a)
{
	for ( Array< real >::iterator it = coordinates_.begin();
	      it != coordinates_.end();
	      ++it )
  {
		*it *= a;
  }

	return *this;
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator/=(real const a)
{
  if(small(a))
  {
    error("MeshGeometry : dividing coordinates by zero");
  }
  real const b = 1.0 / a;

  for ( Array< real >::iterator it = coordinates_.begin();
        it != coordinates_.end();
        ++it )
  {
    *it *= b;
  }
  return *this;
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator+=(real const a)
{
  for ( Array< real >::iterator it = coordinates_.begin();
        it != coordinates_.end();
        ++it )
  {
    *it += a;
  }
  return *this;
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator-=(real const a)
{
  for ( Array< real >::iterator it = coordinates_.begin();
        it != coordinates_.end();
        ++it )
  {
    *it -= a;
  }

  return *this;
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator+=(Point const& p)
{
  for ( Array< real >::iterator it = coordinates_.begin();
        it != coordinates_.end(); )
  {
    for (uint i = 0; i < dim_; ++i)
    {
      *it += p[i];
      ++it;
    }
  }
  return *this;
}
//-----------------------------------------------------------------------------
MeshGeometry& MeshGeometry::operator-=(Point const& p)
{
  for ( Array< real >::iterator it = coordinates_.begin();
        it != coordinates_.end(); )
  {
    for (uint i = 0; i < dim_; ++i)
    {
      *it += p[i];
      ++it;
    }
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
  prm("dimension" , dim_);
  prm("size"      , size_);
  end();
}
//-----------------------------------------------------------------------------
void MeshGeometry::dump() const
{
  for (uint i = 0; i < size_; ++i)
  {
    cout << i << ":";
    for (uint d = 0; d < dim_; ++d)
    {
      cout << " " << x(i)[d];
    }
    cout << "\n";
  }
}
//-----------------------------------------------------------------------------
MeshGeometry const& MeshGeometry::operator>>(Array<real>& A) const
{
  A.assign(coordinates_.data(), coordinates_.data() + dim_ * size_);
  A %= dim_;
  return *this;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

