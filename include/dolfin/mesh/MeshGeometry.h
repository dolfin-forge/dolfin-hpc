// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-08
// Last changed: 2007-11-30

#ifndef __MESH_GEOMETRY_H
#define __MESH_GEOMETRY_H

#include "Point.h"

#include <dolfin/common/types.h>

namespace dolfin
{

/// MeshGeometry stores the geometry imposed on a mesh. Currently,
/// the geometry is represented by the set of coordinates for the
/// vertices of a mesh, but other representations are possible.

class MeshGeometry
{
public:

  /// Create empty set of coordinates
  MeshGeometry();

  /// Copy constructor
  MeshGeometry(MeshGeometry const& geometry);

  /// Destructor
  ~MeshGeometry();

  /// Assignment
  MeshGeometry const& operator=(MeshGeometry const& geometry);

  /// Return Euclidean dimension of coordinate system
  uint dim() const;

  /// Return number of coordinates
  uint size() const;

  /// Return value of coordinate n in direction i
  real& x(uint n, uint i);

  /// Return value of coordinate n in direction i
  real x(uint n, uint i) const;

  /// Return array of values for coordinate n
  real* x(uint n);

  /// Return array of values for coordinate n
  real const * x(uint n) const;

  /// Return array of values for all coordinates
  real* coordinates();

  /// Return array of values for all coordinates
  real const * coordinates() const;

  /// Return coordinate n as a 3D point value
  Point point(uint n) const;

  /// Clear all data
  void clear();

  /// Initialize coordinate list to given geometrical dimension and size
  void init(uint gdim, uint size);

  /// Set value of coordinate n in direction i
  void set(uint n, uint i, real x);

  /// Set value of coordinates n
  void set(uint n, real const * x);

  /// Return token identifying the internal state of mesh geometry
  int token() const;

  /// Display data
  void disp() const;

private:

  // Euclidean dimension
  uint dim_;

  // Number of coordinates
  uint size_;

  // Coordinates for all vertices stored as a contiguous array
  real * coordinates_;

  //
  int timestamp_;

};

//--- INLINES -----------------------------------------------------------------

inline uint MeshGeometry::dim() const
{
  return dim_;
}

//-----------------------------------------------------------------------------
inline uint MeshGeometry::size() const
{
  return size_;
}

//-----------------------------------------------------------------------------
inline real& MeshGeometry::x(uint n, uint i)
{
  dolfin_assert(n < size_ && i < dim_);
  return coordinates_[n * dim_ + i];
}

//-----------------------------------------------------------------------------
inline real MeshGeometry::x(uint n, uint i) const
{
  dolfin_assert(n < size_ && i < dim_);
  return coordinates_[n * dim_ + i];
}

//-----------------------------------------------------------------------------
inline real* MeshGeometry::x(uint n)
{
  return coordinates_ + n * dim_;
}

//-----------------------------------------------------------------------------
inline real const * MeshGeometry::x(uint n) const
{
  return coordinates_ + n * dim_;
}

//-----------------------------------------------------------------------------
inline real* MeshGeometry::coordinates()
{
  return coordinates_;
}

//-----------------------------------------------------------------------------
inline real const * MeshGeometry::coordinates() const
{
  return coordinates_;
}

}

#endif
