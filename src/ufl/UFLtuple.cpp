// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-28
// Last changed: 2014-01-28

#include <dolfin/ufl/UFLtuple.h>

namespace ufl
{

//-----------------------------------------------------------------------------
Object::repr_t const tuple::repr() const
{
  std::stringstream ss;
  ss << "(" << obj_.repr() << ")";
  return ss.str();
}

//-----------------------------------------------------------------------------
std::string const tuple::str() const
{
  std::stringstream ss;
  return "(" + obj_.str() + ")";
}

//-----------------------------------------------------------------------------
void tuple::display() const
{
  std::cout << "Tuple of " << std::endl;
  Object::display();
}

//-----------------------------------------------------------------------------
Object::repr_t const tuple::make_repr(
    std::vector<Object const *> const& prototype) const
{
  std::stringstream ss;
  ss << "(" << Object::make_repr(prototype) << ")";
  return ss.str();
}

} /* namespace icorne */
