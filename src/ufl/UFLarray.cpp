// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-28
// Last changed: 2014-01-28

#include <dolfin/ufl/UFLarray.h>

namespace ufl
{

//-----------------------------------------------------------------------------
array::array(Object const& other, bool const unpack) :
    obj_(other),
    unpack_(unpack)
{
}

//-----------------------------------------------------------------------------
array::~array()
{
}

//-----------------------------------------------------------------------------
Object::repr_t const array::repr() const
{
  std::stringstream ss;
  if (unpack_)
  {
    ss << "*";
  }
  ss << "[" << obj_.repr() << "]";
  return ss.str();
}

//-----------------------------------------------------------------------------
std::string const array::str() const
{
  std::stringstream ss;
  if (unpack_)
  {
    ss << "*";
  }
  ss << "[" << obj_.str() << "]";
  return ss.str();
}

//-----------------------------------------------------------------------------
void array::display() const
{
  std::cout << "Array of ";
  Object::display();
}

//-----------------------------------------------------------------------------
Object::repr_t const array::make_repr(
    std::vector<Object const *> const& prototype) const
{
  std::stringstream ss;
  ss << "[" << Object::make_repr(prototype) << "]";
  return ss.str();
}

//-----------------------------------------------------------------------------
Object::repr_t array::unpack(Object::repr_t const& repr_array)
{
  //FIXME: improve checking of string
  size_t openpos = repr_array.find("[");
  size_t closepos = repr_array.rfind("]");
  return repr_array.substr(openpos + 1, closepos - openpos - 1);
}

} /* namespace icorne */
