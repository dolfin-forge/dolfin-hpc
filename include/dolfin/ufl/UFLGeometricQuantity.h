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
  ~GeometricQuantity();

  ///
  Cell const& cell();

  ///
  virtual ValueArray const& shape() const = 0;

  /// __repr__
  virtual repr_t const& repr() const = 0;

  /// __str__
  virtual std::string const& str() const = 0;

  ///
  virtual void display() const;

private:

  Cell const& cell_;

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_GEOMETRIC_QUANTITY_H */
