// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLClass.h>

#include <iostream>
#include <sstream>

namespace ufl
{

Object::repr_t const Class::default_repr_ = Object::repr_t("Class");
std::string const Class::default_str_ = "class";

//-----------------------------------------------------------------------------
Class::Class() :
    name_()
{
}

//-----------------------------------------------------------------------------
Class::Class(std::string const& name) :
    name_(name)
{
}

//-----------------------------------------------------------------------------
Class::~Class()
{
}

//-----------------------------------------------------------------------------
std::string const& Class::name() const
{
  return name_;
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
Object::repr_t const Class::make_repr(
    std::vector<Object const *> const& prototype) const
{
  std::stringstream ret;
  ret << name() << "(";
  ret << Object::make_repr(prototype);
  ret << ")";
  return ret.str();
}

//-----------------------------------------------------------------------------
void Class::display() const
{
  std::cout << "Class '" << this->name() << "'" << std::endl;
  std::cout << ".str : " << this->str() << std::endl;
  std::cout << ".repr: " << this->repr() << std::endl;
}

//-----------------------------------------------------------------------------
ValueArray::ValueArray() :
    std::vector<uint>()
{
}

//-----------------------------------------------------------------------------
ValueArray::ValueArray(uint const i) :
    std::vector<uint>(1, i)
{
}

//-----------------------------------------------------------------------------
ValueArray::ValueArray(uint const k, uint const i) :
    std::vector<uint>(k, i)
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
  for (ValueArray::const_iterator it = this->begin(); it != this->end(); ++it)
  {
    ss << *it << ",";
  }
  ss << ")";
  return ss.str();
}

}
