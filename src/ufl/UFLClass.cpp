// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLClass.h>

#include <sstream>

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
