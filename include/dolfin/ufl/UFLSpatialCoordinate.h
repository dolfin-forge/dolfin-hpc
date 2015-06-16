// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_SPATIAL_COORDINATE_H_
#define __UFL_SPATIAL_COORDINATE_H_

#include <dolfin/ufl/UFLGeometricQuantity.h>

#include <dolfin/ufl/UFLValueArray.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class
 *
 *  @brief
 */

class Cell;

class SpatialCoordinate : public GeometricQuantity
{

public:

  ///
  SpatialCoordinate(Cell const& cell);

  ///
  ~SpatialCoordinate();

  ///
  bool is_cellwise_constant();

  ///
  ValueArray const& shape() const;

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

private:

  ValueArray const shape_;

  repr_t const repr_;
  std::string const str_;

};

} /* namespace ufl */
#endif /* __UFL_SPATIAL_COORDINATE_H_ */
