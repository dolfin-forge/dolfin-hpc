// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_SPACE_H_
#define __UFL_SPACE_H_

#include <dolfin/ufl/UFLClass.h>

#include <dolfin/common/types.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Cell
 *
 *  @brief  Provides an interface complying with UFL Space.
 */

class Space : public Class
{

public:

  ///
  Space(uint const& dim);

  ///
  ~Space();

  /// Return number of space dimensions
  uint dimension() const;

  /// __repr__
  repr_t const repr() const;

  /// __str__
  std::string const str() const;

private:

  uint const dimension_;

  mutable repr_t repr_;
  mutable std::string str_;

};

} /* namespace ufl */
#endif /* __UFL_SPACE_H */
