// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLClass.h>

namespace dolfin
{

std::string const UFLClass::default_repr_ = "object";
std::string const UFLClass::default_str_ = "Class";

//-----------------------------------------------------------------------------
UFLClass::UFLClass()
{
}

//-----------------------------------------------------------------------------
UFLClass::~UFLClass()
{
}

//-----------------------------------------------------------------------------
std::string const UFLClass::repr() const
{
  return default_repr_;
}

//-----------------------------------------------------------------------------
std::string const UFLClass::str() const
{
  return default_str_;
}

}
