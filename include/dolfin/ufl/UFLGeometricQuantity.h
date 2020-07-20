// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_GEOMETRIC_QUANTITY_H
#define __DOLFIN_UFL_GEOMETRIC_QUANTITY_H

#include <dolfin/ufl/UFLClass.h>

#include <vector>

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
class ValueArray;

class GeometricQuantity : public Class
{
public:

  ///
  GeometricQuantity(std::string const& name, Cell const& cell);

  ///
  ~GeometricQuantity() override;

  ///
  Cell const& cell();

  ///
  virtual ValueArray const& shape() const = 0;

  /// __repr__
  repr_t const& repr() const override = 0;

  /// __str__
  std::string const& str() const override = 0;

  ///
  void display() const override;

private:

  Cell const& cell_;

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_GEOMETRIC_QUANTITY_H */
