// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_GEOMETRY_H
#define __DOLFIN_MESH_GEOMETRY_H

#include <dolfin/common/Tokenized.h>
#include <dolfin/common/Clonable.h>

#include <dolfin/common/types.h>
#include <dolfin/mesh/Point.h>

namespace dolfin
{

template<class T> class Array;
class Mesh;
struct Space;

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

  /// Create a set of point coordinates with given size
  MeshGeometry(Space const& space, uint size = 0);

  /// Copy constructor
  MeshGeometry(MeshGeometry const& other);

  /// Destructor
  ~MeshGeometry();

  /// Assignment
  MeshGeometry & operator=(MeshGeometry const& geometry);

  /// Equality
  bool operator==(MeshGeometry const& other) const;

  /// Non-equality
  bool operator!=(MeshGeometry const& other) const;

  /// Swap instances
  friend void swap( MeshGeometry& a, MeshGeometry& b );

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

  /// Resize space coordinates to size
  void resize(uint size);

  ///
  void finalize();

  /// Set value of all coordinates
  void set(real const * x);

  /// Set absolute geometric tolerance for given topological dimension
  /// The absolute value of the parameter is set as tolerance.
  void set_abs_tolerance(uint dim, real atol);

  /// Set value of coordinates of point n
  void set(uint n, real const * x);

  /// Get value of coordinates of point n
  void get(uint n, real * x) const;

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

  /// Dump data
  void dump() const;

  //--- SERIALIZATION ---------------------------------------------------------
  MeshGeometry const& operator>>(Array<real>& A) const;

  //--- TOKENIZED -------------------------------------------------------------

  /// Return token identifying the internal state of mesh geometry
  int token() const;

private:

  /// Update token value
  void update_token();

  //---------------------------------------------------------------------------

private:

  //
  Space const * space_;

  /// Euclidean dimension
  uint dim_;

  /// Number of coordinates
  uint size_;

  /// Coordinates for all vertices stored as a contiguous array
  Array<real> coordinates_; //!< @todo this should actually be Point
			    //!instead of real

  /// Absolute tolerances
  Array<real> abs_tol_;

  //
  int timestamp_;

};

//--- INLINES -----------------------------------------------------------------

inline real* MeshGeometry::x(uint n)
{
  dolfin_assert(n < size_);
  return coordinates_.data() + n * dim_;
}

//-----------------------------------------------------------------------------
inline real const * MeshGeometry::x(uint n) const
{
  dolfin_assert(n < size_);
  return coordinates_.data() + n * dim_;
}

//-----------------------------------------------------------------------------
inline void MeshGeometry::set(uint n, real const * x)
{
  dolfin_assert(n < size_);
  std::copy(x, x + dim_, coordinates_.data() + n * dim_);
}

//-----------------------------------------------------------------------------
inline void MeshGeometry::get(uint n, real * x) const
{
  dolfin_assert(n < size_);
  real const * xn = coordinates_.data() + n * dim_; std::copy(xn, xn + dim_, x);
}

//-----------------------------------------------------------------------------
inline void MeshGeometry::set(real const * x)
{
  std::copy(x, x + dim_ * size_, coordinates_.data());
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_GEOMETRY_H */
