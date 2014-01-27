// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLClass.h>

#include <sstream>

namespace ufl
{

Object::repr_t const Class::default_repr_ = Object::repr_t("Class");
std::string const Class::default_str_ = "class";

//-----------------------------------------------------------------------------
Class::Class()
{
}

//-----------------------------------------------------------------------------
Class::~Class()
{
}

//-----------------------------------------------------------------------------
Object::repr_t const Class::repr() const
{
  return default_repr_;
}

//-----------------------------------------------------------------------------
std::string const Class::str() const
{
  return default_str_;
}

//-----------------------------------------------------------------------------
ValueArray::ValueArray() : std::vector<uint>()
{
}

//-----------------------------------------------------------------------------
ValueArray::ValueArray(uint const i) : std::vector<uint>(1, i)
{
}

//-----------------------------------------------------------------------------
ValueArray::ValueArray(uint const k, uint const i) : std::vector<uint>(k,i)
{
}

//-----------------------------------------------------------------------------
ValueArray::~ValueArray()
{
}

//-----------------------------------------------------------------------------
std::string const ValueArray::str() const
{
  std::stringstream ss;
  ss << "(";
  for(ValueArray::const_iterator it = this->begin(); it != this->end(); ++it)
  {
    ss << *it << ",";
  }
  ss << ")";
  return ss.str();
}

}
