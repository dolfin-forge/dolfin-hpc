// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-08
// Last changed: 2007-11-30

#ifndef __DOLFIN_MESH_GEOMETRY_H
#define __DOLFIN_MESH_GEOMETRY_H

#include <dolfin/common/Tokenized.h>
#include <dolfin/common/Clonable.h>

#include <dolfin/common/types.h>
#include <dolfin/mesh/Point.h>

namespace dolfin
{

template<class T> class Array;
class Space;

/**
 *  @class  MeshGeometry
 *
 *  @brief  Stores the geometry imposed on a mesh, represented by the set of
 *          coordinates for the vertices of a mesh.
 *
 */

class MeshGeometry : public Tokenized, public Clonable<MeshGeometry>
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

  /// Equality
  bool operator==(MeshGeometry const& other) const;

  /// Non-equality
  bool operator!=(MeshGeometry const& other) const;

  /// Swap instances
  void swap(MeshGeometry& other);

  /// Return space of coordinate system
  Space const& space() const;

  /// Return Euclidean dimension of coordinate system
  uint dim() const;

  /// Return number of points (not coordinates!)
  uint size() const;

  /// Return absolute geometric tolerance for given topological dimension
  real abs_tolerance(uint dim) const;

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

  /// Initialize coordinate list to given geometric dimension and size
  void init(Space const& space, uint size);

  /// Clear all data
  void clear();

  ///
  void finalize();

  /// Set absolute geometric tolerance for given topological dimension
  /// The absolute value of the parameter is set as tolerance.
  void set_abs_tolerance(uint dim, real atol);

  /// Set value of coordinates of point n
  void set(uint n, real const * x);

  /// Assign value of all coordinates
  void assign(real const * x);

  /// Assign values of coordinates from array, size of input should match
  void assign(Array<real> const& coordinates);

  /// Remap coordinates from old to new ordering
  /// The mapping should have the same size as the number of coordinates
  void remap(Array<uint> const& mapping);

  /// Assign coordinates from given mesh geometry and mapping from self to other
  void assign(MeshGeometry const& other, Array<uint> const& mapping);

  /// Scale geometry
  MeshGeometry& operator*=(real const a);
  MeshGeometry& operator/=(real const a);

  /// Offset geometry
  MeshGeometry& operator+=(real const a);
  MeshGeometry& operator-=(real const a);

  /// Translate geometry
  MeshGeometry& operator+=(Point const& p);
  MeshGeometry& operator-=(Point const& p);

  /// Display data
  void disp() const;

  //--- SERIALIZATION ---------------------------------------------------------
  MeshGeometry const& operator>>(Array<real>& A) const;

  //--- TOKENIZED -------------------------------------------------------------

public:

  /// Return token identifying the internal state of mesh geometry
  int token() const;

private:

  /// Update token value
  void update_token();

  //---------------------------------------------------------------------------

private:

  //
  Space const * space_;

  // Euclidean dimension
  uint dim_;

  // Number of coordinates
  uint size_;

  // Coordinates for all vertices stored as a contiguous array
  real * coordinates_;

  // Absolute tolerances
  real * abs_tol_;

  //
  int timestamp_;

};

//--- INLINES -----------------------------------------------------------------

inline real* MeshGeometry::x(uint n)
{
  dolfin_assert(n < size_);
  return coordinates_ + n * dim_;
}

//-----------------------------------------------------------------------------
inline real const * MeshGeometry::x(uint n) const
{
  dolfin_assert(n < size_);
  return coordinates_ + n * dim_;
}

//-----------------------------------------------------------------------------
inline void MeshGeometry::set(uint n, real const * x)
{
  dolfin_assert(n < size_);
  std::copy(x, x + dim_, coordinates_ + n * dim_);
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_GEOMETRY_H */
