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
  ~Circumradius();

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
#endif /* __DOLFIN_UFL_CIRCUMRADIUS_H */
