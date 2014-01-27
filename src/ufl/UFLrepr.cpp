// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-27
// Last changed: 2014-01-27

#include <dolfin/ufl/UFLrepr.h>

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLObject.h>

#include <sstream>

namespace ufl
{

//-----------------------------------------------------------------------------
repr::repr() :
    std::string()
{
}

//-----------------------------------------------------------------------------
repr::repr(std::string const& s) :
    std::string(s)
{
}

//-----------------------------------------------------------------------------
repr::repr(Class const& owner, Object const& arg1) :
    std::string(make_representation(owner, arg1))
{
}

//-----------------------------------------------------------------------------
repr::repr(Class const& owner, Object const& arg1, Object const& arg2) :
    std::string(make_representation(owner, arg1, arg2))
{
}

//-----------------------------------------------------------------------------
repr::repr(Class const& owner, std::vector<Object const *> const& prototype) :
    std::string(make_representation(owner, prototype))
{
}

//-----------------------------------------------------------------------------
repr::~repr()
{
}

//-----------------------------------------------------------------------------
std::string const repr::make_representation(Class const& owner,
                                            Object const& arg1)
{
  std::stringstream ret;
  ret << owner.name() << "(" << arg1.repr() << ")";
  return ret.str();
}

//-----------------------------------------------------------------------------
std::string const repr::make_representation(Class const& owner,
                                            Object const& arg1,
                                            Object const& arg2)
{
  std::stringstream ret;
  ret << owner.name() << "(" << arg1.repr() << ", " << arg2.repr() << ")";
  return ret.str();
}

//-----------------------------------------------------------------------------
std::string const repr::make_representation(
    Class const& owner, std::vector<Object const *> const& prototype)
{
  std::stringstream ret;
  ret << owner.name() << "(";
  std::vector<Object const *>::const_iterator arg = prototype.begin();
  for (; arg != prototype.end(); ++arg, ret << ", ")
  {
    ret << (*arg)->repr();
  }
  ret << ")";
  return ret.str();
}

} /* namespace icorne */
