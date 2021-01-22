// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.

// Original documentation and comments were incorrect, class rewritten in 2015.

#ifndef __DOLFIN_EQUI_AFFINE_MAPPING_H
#define __DOLFIN_EQUI_AFFINE_MAPPING_H

#include <dolfin/mesh/mappings/Mapping.h>

#include <dolfin/mesh/Point.h>

namespace dolfin
{

class Cell;
class Mesh;

/// This class represents the affine map from regular simplices to the current
/// element.
///
/// In 2D, the equilateral triangle is given by:
///
///   {
///     ( 0         , 0   ),
///     ( sqrt(3)/2 , 0.5 ),
///     ( 0         , 1   )
///   }
///
/// In 3D, the regular tetrahedron is given by:
///
///   {
///     ( 0         , 0   , 0         ),
///     ( sqrt(3)/2 , 0.5 , 0         ),
///     ( 0         , 1   , 0         ),
///     ( sqrt(3)/6 , 0.5 , sqrt(2/3) )
///   }
///
/// The dimension d of the map is automatically determined from the
/// arguments used when calling the map.

class EquiAffineMapping : public Mapping
{

public:
  /// Constructor
  EquiAffineMapping( Mesh const & mesh );

  /// Destructor
  ~EquiAffineMapping() override;

  /// Update map for current element
  void update( Cell const & cell ) override;

  /// Map given point from the reference element
  void map_from_reference_cell( real const * xref, real * x ) const override;

  /// Map given point to the reference element
  void map_to_reference_cell( real const * x, real * xref ) const override;

  // Determinant of Jacobian of map
  real det;

  // Jacobian of map
  real * J;

  // Inverse of Jacobian of map
  real * K;

private:
  // Update affine map from reference triangle
  void updateTriangle( Cell const & cell );

  // Update affine map from reference tetrahedron
  void updateTetrahedron( Cell const & cell );

  //
  static size_t const d_ = Point::MAX_SIZE;
  size_t const        gdim_;

  // Vertices of current cell
  real p0[d_];
  real p1[d_];
  real p2[d_];
  real p3[d_];
};

}

#endif
