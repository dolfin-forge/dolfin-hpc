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
repr::repr(const char*& s) :
    std::string(s)
{
}

//-----------------------------------------------------------------------------
repr::repr(Class const& owner, Object const& arg1) :
    std::string(make_repr(owner, arg1))
{
}

//-----------------------------------------------------------------------------
repr::repr(Class const& owner, Object const& arg1, Object const& arg2) :
    std::string(make_repr(owner, arg1, arg2))
{
}

//-----------------------------------------------------------------------------
repr::repr(Class const& owner, std::vector<Object const *> const& prototype) :
    std::string(owner.make_repr(prototype))
{
}

//-----------------------------------------------------------------------------
repr::~repr()
{
}

//-----------------------------------------------------------------------------
std::string const repr::make_repr(Class const& owner, Object const& arg1)
{
  std::vector<Object const *> p;
  p.push_back(&arg1);
  return owner.make_repr(p);
}

//-----------------------------------------------------------------------------
std::string const repr::make_repr(Class const& owner, Object const& arg1,
                                  Object const& arg2)
{
  std::vector<Object const *> p;
  p.push_back(&arg1);
  p.push_back(&arg2);
  return owner.make_repr(p);
}

} /* namespace icorne */
