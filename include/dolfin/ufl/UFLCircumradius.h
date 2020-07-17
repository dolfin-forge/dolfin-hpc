// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_CIRCUMRADIUS_H
#define __DOLFIN_UFL_CIRCUMRADIUS_H

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

class Circumradius : public GeometricQuantity
{

public:

  ///
  Circumradius(Cell const& cell);

  ///
  ~Circumradius() override;

  ///
  ValueArray const& shape() const override;

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

private:

  ValueArray const shape_;

  repr_t const repr_;
  std::string const str_;

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_CIRCUMRADIUS_H */
