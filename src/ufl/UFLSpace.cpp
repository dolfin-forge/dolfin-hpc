// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLSpace.h>

#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLSpace::UFLSpace(uint const& dim) :
    dimension_(dim)
{
  std::stringstream ssrepr;
  ssrepr << "Space(" << dimension_ << ")";
  repr_ = ssrepr.str();

  std::stringstream ssstr;
  ssstr << "R" << dimension_;
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
UFLSpace::~UFLSpace()
{
}

//-----------------------------------------------------------------------------
uint UFLSpace::dimension() const
{
  return dimension_;
}

//-----------------------------------------------------------------------------
std::string const UFLSpace::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const UFLSpace::str() const
{
  return str_;
}

}
