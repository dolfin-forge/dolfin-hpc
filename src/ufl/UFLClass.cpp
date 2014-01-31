// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLClass.h>

#include <dolfin/log/log.h>

#include <iostream>
#include <sstream>

namespace ufl
{

Object::repr_t const Class::default_repr_ = Object::repr_t("Class");
std::string const Class::default_str_ = "class";

//-----------------------------------------------------------------------------
Class::Class() :
    args_repr_(),
    name_()
{
}

//-----------------------------------------------------------------------------
Class::Class(std::string const& name) :
    cpp_proto_(name, std::vector<Object const *>()),
    args_repr_(),
    name_(name)
{
}

//-----------------------------------------------------------------------------
Class::Class(std::string const& name, repr_t const& repr) :
    cpp_proto_(make_proto(repr)),
    args_repr_(make_args_repr(repr)),
    name_(name)
{
  if (name_ != cpp_proto_.first)
  {
    dolfin::error(
        "The representation is invalid as it does not match the class");
  }
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
    std::vector<Object const *> const& args) const
{
  std::stringstream ret;
  ret << name() << "(";
  ret << Object::make_repr(args);
  ret << ")";
  return ret.str();
}

//-----------------------------------------------------------------------------
Object::repr_t const Class::make_repr(Object const *& arg1) const
{
  std::vector<Object const *> args;
  args.push_back(arg1);
  return Object::make_repr(args);
}

//-----------------------------------------------------------------------------
Object::repr_t const Class::make_repr(Object const *& arg1,
                                      Object const *& arg2) const
{
  std::vector<Object const *> args;
  args.push_back(arg1);
  args.push_back(arg2);
  return Object::make_repr(args);
}

//-----------------------------------------------------------------------------
std::vector<Object::repr_t> const Class::make_args_repr(
    repr_t const& repr) const
{
  std::string str = repr;
  std::string dopen = "(";
  std::string dclose = ")";

  size_t openpos = str.find(dopen);
  size_t closepos = str.rfind(dclose);
  std::string name = str.substr(0, openpos);
//  std::cout << openpos << "-" << closepos << std::endl;
//  std::cout << "arg = " << str.substr(openpos + 1, closepos - openpos - 1)
//      << std::endl;
  std::vector<repr_t> args = Object::make_args_repr(
      str.substr(openpos + 1, closepos - openpos - 1));
  return args;
}

//-----------------------------------------------------------------------------
Class::CppProto Class::make_proto(Object::repr_t repr) const
{
  std::string str = repr;
  size_t openpos = str.find("(");
  std::string name = str.substr(0, openpos);
  str.erase(0, openpos);
  std::vector<repr_t> args = Class::make_args_repr(str);
  return CppProto(name, Object::make_args(args));
}

//-----------------------------------------------------------------------------
Object::repr_t const& Class::arg(size_t i)
{
  if(i >= args_repr_.size())
  {
    for(size_t k=0; k < args_repr_.size(); ++k)
    {
      std::cout << args_repr_[k] << ", ";
    }
    dolfin::error("Invalid index required for argument.");
  }
  return args_repr_[i];
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
