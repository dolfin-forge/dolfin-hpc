// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __UFL_EQUATION_H_
#define __UFL_EQUATION_H_

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
  ~Equation();

  //--- INTERFACE inherited from UFLClass -------------------------------------

  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:
  Form const lhs_;
  Form const rhs_;

  repr_t const repr_;
  std::string const str_;
};

} /* namespace ufl */
#endif /* __UFL_EQUATION_H_ */
