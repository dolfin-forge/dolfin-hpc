// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLQuadratureScheme.h>

#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLQuadratureScheme::UFLQuadratureScheme() :
    UFLClass(),
    repr_("QuadratureScheme"),
    str_()
{
}

//-----------------------------------------------------------------------------
UFLQuadratureScheme::~UFLQuadratureScheme()
{
}

//-----------------------------------------------------------------------------
std::string const UFLQuadratureScheme::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const UFLQuadratureScheme::str() const
{
  return str_;
}

}
