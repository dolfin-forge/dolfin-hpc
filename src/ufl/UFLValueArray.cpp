// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLValueArray.h>
#include <sstream>

namespace ufl
{

//-----------------------------------------------------------------------------
ValueArray::ValueArray() :
    std::vector<dolfin::uint>()
{
}

//-----------------------------------------------------------------------------
ValueArray::ValueArray(dolfin::uint const i) :
    std::vector<dolfin::uint>(1, i)
{
}

//-----------------------------------------------------------------------------
ValueArray::ValueArray(dolfin::uint const k, dolfin::uint const i) :
    std::vector<dolfin::uint>(k, i)
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

//-----------------------------------------------------------------------------
dolfin::uint ValueArray::prod() const
{
  if(this->size() == 0)
  {
    return 0;
  }
  dolfin::uint vs = 1;
  for (ValueArray::const_iterator j = this->begin(); j != this->end(); ++j)
  {
    vs *= *j;
  }
  return vs;
}

//-----------------------------------------------------------------------------
ValueArray operator+(ValueArray const&v1, ValueArray const& v2)
{
  ValueArray ret(v1);
  ret.reserve(v1.size() + v2.size());
  ret.insert(ret.end(), v2.begin(), v2.end());
  return ret;
}

}
