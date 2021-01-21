// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_GEOMETRY_H
#define __DOLFIN_MESH_GEOMETRY_H

#include <dolfin/common/Array.h>
#include <dolfin/common/Clonable.h>
#include <dolfin/common/Tokenized.h>
#include <dolfin/common/types.h>
#include <dolfin/mesh/Point.h>

namespace dolfin
{

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
  ~MeshGeometry() override;

  /// Assignment
  auto operator=(MeshGeometry const& geometry) -> MeshGeometry &;

  /// Equality
  auto operator==(MeshGeometry const& other) const -> bool;

  /// Non-equality
  auto operator!=(MeshGeometry const& other) const -> bool;

  /// Swap instances
  friend void swap( MeshGeometry& a, MeshGeometry& b );

  /// Return space of coordinate system
  auto space() const -> Space const&;

  /// Return Euclidean dimension of coordinate system
  auto dim() const -> uint;

  /// Return number of points (not coordinates!)
  auto size() const -> uint;

  /// Return absolute geometric tolerance for given topological dimension
  auto abs_tolerance(uint dim) const -> real;

  /// Return array of values for coordinate n
  auto x(uint n) -> real*;

  /// Return array of values for coordinate n
  auto x(uint n) const -> real const *;

  /// Return array of values for all coordinates
  auto coordinates() -> real*;

  /// Return array of values for all coordinates
  auto coordinates() const -> real const *;

  /// Return coordinate n as a 3D point value
  auto point(uint n) const -> Point;

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
  auto operator*=(real const a) -> MeshGeometry&;
  auto operator/=(real const a) -> MeshGeometry&;

  /// Offset geometry
  auto operator+=(real const a) -> MeshGeometry&;
  auto operator-=(real const a) -> MeshGeometry&;

  /// Translate geometry
  auto operator+=(Point const& p) -> MeshGeometry&;
  auto operator-=(Point const& p) -> MeshGeometry&;

  /// Display data
  void disp() const;

  /// Dump data
  void dump() const;

  //--- SERIALIZATION ---------------------------------------------------------
  auto operator>>(Array<real>& A) const -> MeshGeometry const&;

  //--- TOKENIZED -------------------------------------------------------------

  /// Return token identifying the internal state of mesh geometry
  auto token() const -> int override;

private:

  /// Update token value
  void update_token() override;

  //---------------------------------------------------------------------------

private:

  //
  Space const * space_;

  /// Euclidean dimension
  uint dim_;

  /// Number of coordinates
  uint size_;

  /// Coordinates for all vertices stored as a contiguous array
  Array<real> coordinates_; //!< @todo this should actually be Point
			    //!instead of real

  /// Absolute tolerances
  Array<real> abs_tol_;

  //
  int timestamp_;

};

//--- INLINES -----------------------------------------------------------------

inline auto MeshGeometry::x( uint n ) -> real *
{
  dolfin_assert( n < size_ );
  return coordinates_.data() + n * dim_;
}

//-----------------------------------------------------------------------------
inline auto MeshGeometry::x( uint n ) const -> real const *
{
  dolfin_assert( n < size_ );
  return coordinates_.data() + n * dim_;
}

//-----------------------------------------------------------------------------
inline void MeshGeometry::set( uint n, real const * x )
{
  dolfin_assert( n < size_ );
  std::copy( x, x + dim_, coordinates_.data() + n * dim_ );
}

//-----------------------------------------------------------------------------
inline void MeshGeometry::get( uint n, real * x ) const
{
  dolfin_assert( n < size_ );
  real const * xn = coordinates_.data() + n * dim_;
  std::copy( xn, xn + dim_, x );
}

//-----------------------------------------------------------------------------
inline void MeshGeometry::set( real const * x )
{
  std::copy( x, x + dim_ * size_, coordinates_.data() );
}

//-----------------------------------------------------------------------------
inline auto MeshGeometry::space() const -> Space const &
{
  return *space_;
}

//-----------------------------------------------------------------------------
inline auto MeshGeometry::dim() const -> uint
{
  return dim_;
}

//-----------------------------------------------------------------------------
inline auto MeshGeometry::size() const -> uint
{
  return size_;
}

//-----------------------------------------------------------------------------
inline auto MeshGeometry::abs_tolerance( uint dim ) const -> real
{
  dolfin_assert( dim <= dim_ );
  return abs_tol_[dim];
}

//-----------------------------------------------------------------------------
inline auto MeshGeometry::point( uint n ) const -> Point
{
  return Point( coordinates_.data() + n * dim_ );
}

//-----------------------------------------------------------------------------
inline auto MeshGeometry::coordinates() -> real *
{
  return coordinates_.data();
}

//-----------------------------------------------------------------------------
inline auto MeshGeometry::coordinates() const -> real const *
{
  return coordinates_.data();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_GEOMETRY_H */
