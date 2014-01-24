// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_CIRCUMRADIUS_H_
#define __UFL_CIRCUMRADIUS_H_

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

class UFLCircumradius : public UFLGeometricQuantity
{

public:

  ///
  UFLCircumradius(UFLCell const& cell);

  ///
  ~UFLCircumradius();

  ///
  std::vector<uint> const& shape() const;

  /// __repr__
  std::string const repr() const;

  /// __str__
  std::string const str() const;

private:

  std::vector<uint> const shape_;

  std::string const repr_;
  std::string const str_;

};

} /* namespace dolfin */
#endif /* __UFL_CIRCUMRADIUS_H_ */
