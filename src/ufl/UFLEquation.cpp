// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLEquation.h>

namespace ufl
{

//-----------------------------------------------------------------------------
Equation::Equation(Form const& lhs, Form const& rhs) :
    Class("Equation"),
    lhs_(lhs),
    rhs_(rhs),
    repr_(*this, lhs_, rhs_),
    str_("")
{
}

//-----------------------------------------------------------------------------
Equation::Equation(repr_t const & repr) :
    Class("Equation", repr),
    lhs_(arg(0)),
    rhs_(arg(1)),
    repr_(*this, lhs_, rhs_),
    str_("")
{
}

//-----------------------------------------------------------------------------
Equation::~Equation()
{
}

//-----------------------------------------------------------------------------
Object::repr_t const& Equation::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& Equation::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
void Equation::display() const
{
}

}
