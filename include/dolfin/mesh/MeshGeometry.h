// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-08
// Last changed: 2007-11-30

#ifndef __DOLFIN_MESH_GEOMETRY_H
#define __DOLFIN_MESH_GEOMETRY_H

#include <dolfin/common/Tokenized.h>

#include <dolfin/common/types.h>
#include <dolfin/mesh/EuclideanSpace.h>
#include <dolfin/mesh/Point.h>

namespace dolfin
{

template<class T> class Array;

/// MeshGeometry stores the geometry imposed on a mesh. Currently,
/// the geometry is represented by the set of coordinates for the
/// vertices of a mesh, but other representations are possible.

class MeshGeometry : public Tokenized
{
public:

  /// Create empty set of coordinates
  MeshGeometry();

  /// Create set of coordinates given geometric dimension and size
  MeshGeometry(uint gdim, uint size);

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

  /// Return absolute geometric tolerance for given topological dimension
  real abs_tolerance(uint dim) const;

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

  ///
  void finalize();

  /// Initialize coordinate list to given geometrical dimension and size
  void init(uint gdim, uint size);

  /// Set absolute geometric tolerance for given topological dimension
  /// The absolute value of the parameter is set as tolerance.
  void set_abs_tolerance(uint dim, real atol);

  /// Set value of coordinate n in direction i
  void set(uint n, uint i, real x);

  /// Set value of coordinates n
  void set(uint n, real const * x);

  /// Remap coordinates from old to new ordering
  /// The mapping should have the same size as the number of coordinates
  void remap(Array<uint> const& map);

  /// Display data
  void disp() const;

  //--- CHECK ROUTINES --------------------------------------------------------

  /// Check
  void check() const;

  //--- TOKENIZED -------------------------------------------------------------

public:

  /// Return token identifying the internal state of mesh geometry
  int token() const;

private:

  /// Update token value
  void update_token();

  //---------------------------------------------------------------------------

private:

  // Euclidean dimension
  uint dim_;

  // Number of coordinates
  uint size_;

  // Absolute tolerances
  real abs_tol_[EuclideanSpace::MAX_DIMENSION+1];

  // Coordinates for all vertices stored as a contiguous array
  real * coordinates_;

  //
  int timestamp_;

};

//--- INLINES -----------------------------------------------------------------

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

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_GEOMETRY_H */
