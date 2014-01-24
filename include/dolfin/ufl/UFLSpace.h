// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_SPACE_H_
#define __UFL_SPACE_H_

#include <dolfin/ufl/UFLClass.h>

#include <dolfin/common/types.h>

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLCell
 *
 *  @brief  Provides an interface complying with UFL Space.
 */

class UFLSpace : public UFLClass
{

public:

  ///
  UFLSpace(uint const& dim);

  ///
  ~UFLSpace();

  /// Return number of space dimensions
  uint dimension() const;

  /// __repr__
  std::string const repr() const;

  /// __str__
  std::string const str() const;

private:

  uint const dimension_;

  mutable std::string repr_;
  mutable std::string str_;

};

} /* namespace dolfin */
#endif /* __UFL_SPACE_H */
