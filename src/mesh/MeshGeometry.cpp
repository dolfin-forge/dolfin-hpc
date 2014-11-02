// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-19
// Last changed: 2006-10-19

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/MeshGeometry.h>

#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry() :
    dim_(0),
    size_(0),
    coordinates_(0),
    timestamp_(time(0))
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry(MeshGeometry const& geometry) :
    dim_(0),
    size_(0),
    coordinates_(0),
    timestamp_(time(0))
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
  const uint n = dim_ * size_;
  coordinates_ = new real[n];

  // Copy data
  for (uint i = 0; i < n; ++i)
  {
    coordinates_[i] = geometry.coordinates_[i];
  }

  timestamp_ = geometry.timestamp_;

  return *this;
}
//-----------------------------------------------------------------------------
Point MeshGeometry::point(uint n) const
{
  real _x = 0.0;
  real _y = 0.0;
  real _z = 0.0;

  if (dim_ > 0) _x = x(n, 0);
  if (dim_ > 1) _y = x(n, 1);
  if (dim_ > 2) _z = x(n, 2);

  Point p(_x, _y, _z);
  return p;
}
//-----------------------------------------------------------------------------
void MeshGeometry::clear()
{
  dim_ = 0;
  size_ = 0;
  delete[] coordinates_;
  coordinates_ = NULL;
}
//-----------------------------------------------------------------------------
void MeshGeometry::init(uint gdim, uint size)
{
  // Delete old data if any
  clear();

  // Allocate new data
  coordinates_ = new real[gdim * size];

  // Save dimension and size
  dim_ = gdim;
  size_ = size;
  timestamp_ = time(0);
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
int MeshGeometry::token() const
{
  return timestamp_ + size(); // FIXME
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

}

