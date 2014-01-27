// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_CELL_SURFACE_AREA_H_
#define __UFL_CELL_SURFACE_AREA_H_

#include <dolfin/ufl/UFLGeometricQuantity.h>

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class
 *
 *  @brief
 */

class UFLCell;

class UFLCellSurfaceArea : public UFLGeometricQuantity
{

public:

  ///
  UFLCellSurfaceArea(UFLCell const& cell);

  ///
  ~UFLCellSurfaceArea();

  ///
  ValueArray const& shape() const;

  /// __repr__
  std::string const repr() const;

  /// __str__
  std::string const str() const;

private:

  ValueArray const shape_;

  std::string const repr_;
  std::string const str_;

};

} /* namespace dolfin */
#endif /* __UFL_CELL_SURFACE_AREA_H_ */
