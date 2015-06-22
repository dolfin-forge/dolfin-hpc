// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLQuadratureScheme.h>

#include <sstream>

namespace ufl
{

//-----------------------------------------------------------------------------
QuadratureScheme::QuadratureScheme() :
    Class("QuadratureScheme"),
    repr_("None"),
    str_("None")
{
}

//-----------------------------------------------------------------------------
QuadratureScheme::~QuadratureScheme()
{
}

//-----------------------------------------------------------------------------
Object::repr_t const& QuadratureScheme::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& QuadratureScheme::str() const
{
  return str_;
}

}
