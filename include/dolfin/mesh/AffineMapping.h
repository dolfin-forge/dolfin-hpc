// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU GPL Version 2.
//
// This version fixes several issues compared to the original:
// - assumption that topological dimension is equal to geometric dimension.
// - duplication of code.
//
// First added:
// Last changed:

#ifndef __DOLFIN_AFFINE_MAPPING_H
#define __DOLFIN_AFFINE_MAPPING_H

#include <dolfin/mesh/Mapping.h>

#include <dolfin/common/constants.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/Cell.h>

namespace dolfin
{

/// This class represents the affine map from the reference element to
/// the current element.
///
/// The 2D reference element is given by (0,0)-(1,0)-(0,1).
/// The 3D reference element is given by (0,0,0)-(1,0,0)-(0,1,0)-(0,0,1).
///
/// The dimension d of the map is automatically determined from the
/// arguments used when calling the map.

class AffineMapping : public Mapping
{

public:

  /// Constructor
  AffineMapping(Mesh const& mesh);

  /// Destructor
  ~AffineMapping();

  /// Update map for current element
  void update(Cell const& cell);

  /// Map given point from the reference element
  void map_from_reference_cell(real const * xref, real * x) const;

  /// Map given point to the reference element
  void map_to_reference_cell(real const * x, real * xref) const;

  // Determinant of Jacobian of map
  real det;

  // Jacobian of map
  real *J;

  // Inverse of Jacobian of map
  real *K;

private:

  // Update affine map from reference interval
  void updateInterval(Cell const& cell);

  // Update affine map from reference triangle
  void updateTriangle(Cell const& cell);

  // Update affine map from reference tetrahedron
  void updateTetrahedron(Cell const& cell);

  // Update
  void updateR1();
  void updateR2();
  void updateR3();

  //
  static uint const n_ = 8; // Maximum number of vertices per cell
  static uint const d_ = Point::MAX_SIZE;
  uint const gdim_;

  // Vertices of current cell
  real p[n_][d_];

};

}

#endif
