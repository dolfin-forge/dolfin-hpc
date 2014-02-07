// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLTensors.h>

//#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  ComponentTensor::ComponentTensor(Expression const& expression, IndexBase const& index) :
    Class("ComponentTensor"),
    expression_(expression),
    index_(index)
  {
    std::stringstream ssrepr;
    ssrepr << "ComponentTensor(" << expression_.repr() << "," << index_.repr() << ")";
    repr_ = ssrepr.str();

    //Is this the same implementation as in python?
    std::stringstream ssstr;
    ssstr << "{ A | A_{" << index_.str() << "} = " << expression_.str() << " }";
    str_ = ssstr.str();
  }

//-----------------------------------------------------------------------------
  ComponentTensor::ComponentTensor(repr_t const & repr) :
    Class("ComponentTensor", repr),
    expression_(arg(0)),
    index_(arg(1))
  {
  }

//-----------------------------------------------------------------------------
  ComponentTensor::~ComponentTensor()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression const, MultiIndex const> const& ComponentTensor::operands() const
  {
    return std::make_pair(expression_, index_);  
  }

//-----------------------------------------------------------------------------
  Object::repr_t const ComponentTensor::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const ComponentTensor::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void ComponentTensor::display() const
  {
  }

}
