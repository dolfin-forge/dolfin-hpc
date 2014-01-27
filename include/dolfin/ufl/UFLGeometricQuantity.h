// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_GEOMETRIC_QUANTITY_H_
#define __UFL_GEOMETRIC_QUANTITY_H_

#include <dolfin/ufl/UFLClass.h>

#include <vector>

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

class UFLGeometricQuantity : public UFLClass
{
public:

  ///
  UFLGeometricQuantity(UFLCell const& cell);

  ///
  ~UFLGeometricQuantity();

  ///
  UFLCell const& cell();

  ///
  virtual ValueArray const& shape() const = 0;

  /// __repr__
  virtual std::string const repr() const = 0;

  /// __str__
  virtual std::string const str() const = 0;

private:

  UFLCell const& cell_;

};

} /* namespace dolfin */
#endif /* __UFL_GEOMETRIC_QUANTITY_H_ */
