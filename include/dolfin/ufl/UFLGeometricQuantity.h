// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_GEOMETRIC_QUANTITY_H_
#define __UFL_GEOMETRIC_QUANTITY_H_

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
#endif /* __UFL_GEOMETRIC_QUANTITY_H_ */
