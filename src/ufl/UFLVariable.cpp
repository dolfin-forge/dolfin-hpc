// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLVariable.h>
//#include <dolfin/ufl/UFLAlgebra.h>
//#include <dolfin/ufl/UFLDifferentiation.h>
//#include <dolfin/ufl/UFLExpression.h>
//#include <dolfin/ufl/UFLIndexed.h>
//#include <dolfin/ufl/UFLIndexSum.h>
//#include <dolfin/ufl/UFLTensors.h>

//#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  Label::Label(dolfin::uint const& count) :
    count_(count),
    repr_(*this, count_),
    str_("Label(" + count_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Label::Label(repr_t const& repr) :
    count_(arg(0)),
    repr_(*this, count_),
    str_("Label(" + count_.str() + ")")
  {
  }
//-----------------------------------------------------------------------------
  Label::~Label()
  {
  }
  
//-----------------------------------------------------------------------------
  type<dolfin::uint> const& Label::count() const
  {
    return count_; 
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Label::repr() const
  {
    return repr_;
  }

  //-----------------------------------------------------------------------------
  std::string const Label::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Label::display() const
  {
  }

//-----------------------------------------------------------------------------
  Variable::Variable(Expression const& expression, Label const& label) :
    expression_(expression),
    label_(label),
    repr_(*this, expression_, label_),
    str_("var" + label_.count().str() + "(" + expression_.str() + ")")
  {
  }

//-----------------------------------------------------------------------------
  Variable::Variable(repr_t const& repr) :
    expression_(arg(0)),
    label_(arg(1)),
    repr_(*this, expression_, label_),
    str_("var" + label_.count().str() + "(" + expression_.str() + ")")
  {
  }
//-----------------------------------------------------------------------------
  Variable::~Variable()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression const, Label const> const& Variable::operands() const
  {
    return std::make_pair(expression_, label_);  
  }

//-----------------------------------------------------------------------------
//  Cell const& Variable::cell() const
//  {
//    return expression_.cell();
//  }

//-----------------------------------------------------------------------------
  bool const Variable::is_cellwise_constant() const
  {
    return expression_.is_cellwise_constant();
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Variable::repr() const
  {
    return repr_;
  }

  //-----------------------------------------------------------------------------
  std::string const Variable::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Variable::display() const
  {
  }
}
