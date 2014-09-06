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
    expr_index_(expression, index)
  {
    std::stringstream ssrepr;
    ssrepr << "ComponentTensor(" << expression.repr() << "," << index.repr() << ")";
    repr_ = ssrepr.str();

    //Is this the same implementation as in python?
    std::stringstream ssstr;
    ssstr << "{ A | A_{" << index.str() << "} = " << expression.str() << " }";
    str_ = ssstr.str();
  }

//-----------------------------------------------------------------------------
  ComponentTensor::ComponentTensor(repr_t const & repr) :
    Class("ComponentTensor", repr),
    expr_index_(arg(0), arg(1))
  {
  }

//-----------------------------------------------------------------------------
  ComponentTensor::~ComponentTensor()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression, MultiIndex> const& ComponentTensor::operands() const
  {
    return expr_index_;
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
