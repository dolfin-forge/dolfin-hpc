// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLArgument.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  Argument::Argument(FiniteElementBase const& fe, dolfin::uint const& c) :
    Class("Argument"),
    finite_element_(fe),
    count_(c),
    repr_(*this, finite_element_, count_),
    str_((count_ < 10 ? "v_" + count_.str() : "v_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  Argument::Argument(repr_t const& repr) :
    Class("Argument"),
    finite_element_(*FiniteElementBase::create(arg(0))),
    count_(arg(1)),
    repr_(*this, finite_element_, count_),
    str_((count_ < 10 ? "v_" + count_.str() : "v_{" + count_.str() + "}"))
  {
  }

//-----------------------------------------------------------------------------
  Argument::~Argument()
  {
  }
  
//-----------------------------------------------------------------------------
  FiniteElementBase const& Argument::element() const
  {
    return finite_element_;  
  }

//-----------------------------------------------------------------------------
  ValueArray const& Argument::shape() const
  {
    return finite_element_.value_shape();  
  }

//-----------------------------------------------------------------------------
  Cell const& Argument::cell() const
  {
    return finite_element_.cell();  
  }

//-----------------------------------------------------------------------------
  bool const Argument::is_cellwise_constant() const
  {
    return false;
  }
 
//-----------------------------------------------------------------------------
  Object::repr_t const Argument::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Argument::str() const
  {
    return str_;
  }
 
//-----------------------------------------------------------------------------
  void Argument::display() const
  {
  }

}
