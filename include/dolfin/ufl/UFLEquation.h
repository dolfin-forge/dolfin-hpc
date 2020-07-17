// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_EQUATION_H
#define __DOLFIN_UFL_EQUATION_H

#include <dolfin/ufl/UFLForm.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Equation
 *
 *  @brief  Provides an interface complying with UFL Equation.
 *  This class is used to represent equations expressed by the "=="
 *  operator. Examples include a == L and F == 0 where a, L and F are
 *  Form objects.
 */

class Equation : public Class
{
public:

  ///
  Equation(Form const& lhs, Form const& rhs);

  ///
  Equation(repr_t const& repr);

  ///
  ~Equation() override;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

  ///
  void display() const override;

private:
  Form const lhs_;
  Form const rhs_;

  repr_t const repr_;
  std::string const str_;
};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_EQUATION_H */
