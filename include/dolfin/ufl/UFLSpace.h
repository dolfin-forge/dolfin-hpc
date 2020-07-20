// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_SPACE_H
#define __DOLFIN_UFL_SPACE_H

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLtype.h>

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
  Space(dolfin::uint const& dim);

  ///
  Space(repr_t const& repr);

  ///
  ~Space() override;

  /// UFL: Return number of space dimensions
  dolfin::uint dimension() const;

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

  ///
  void display() const override;

private:

  type<dolfin::uint> const dimension_;

  repr_t const repr_;
  mutable std::string str_;

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_SPACE_H */
