// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_FACET_NORMAL_H_
#define __UFL_FACET_NORMAL_H_

#include <dolfin/ufl/UFLGeometricQuantity.h>

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

class FacetNormal : public GeometricQuantity
{

public:

  ///
  FacetNormal(Cell const& cell);

  ///
  ~FacetNormal();

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

} /* namespace ufl */
#endif /* __UFL_FACET_NORMAL_H_ */
