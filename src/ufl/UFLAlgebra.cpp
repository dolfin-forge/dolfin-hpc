// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLAlgebra.h>

//#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  Sum::Sum(Expression const& s1, Expression const& s2) :
    Class("Sum"),
    s_(s1, s2),
    repr_(*this, s_.first, s_.second),
    str_(s_.first.str() + " + " + s_.second.str())
  {
  }

//-----------------------------------------------------------------------------
  Sum::Sum(repr_t const & repr):
    Class("Sum", repr),
    s_(arg(0),arg(1)),
    repr_(*this, s_.first, s_.second),
    str_(s_.first.str() + " + " + s_.second.str())
  {
  }

//-----------------------------------------------------------------------------
  Sum::~Sum()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression, Expression> const& Sum::operands() const
  {
    return s_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Sum::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Sum::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Sum::display() const
  {
  }

//-----------------------------------------------------------------------------
  Product::Product(Expression const& p1, Expression const& p2) :
    Class("Product"),
    p_(p1, p2),
    repr_(*this, p_.first, p_.second),
    str_(p_.first.str() + " * " + p_.second.str())
  {
  }

//-----------------------------------------------------------------------------
  Product::Product(repr_t const & repr):
    Class("Product", repr),
    p_(arg(0), arg(1)),
    repr_(*this, p_.first, p_.second),
    str_(p_.first.str() + " * " + p_.second.str())
  {
  }

//-----------------------------------------------------------------------------
  Product::~Product()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression, Expression> const& Product::operands() const
  {
    return p_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Product::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Product::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Product::display() const
  {
  }

//-----------------------------------------------------------------------------
  Division::Division(Expression const& d1, Expression const& d2) :
    Class("Division"),
    d_(d1, d2),
    repr_(*this, d_.first, d_.second),
    str_(d_.first.str() + " / " + d_.second.str())
  {
  }

//-----------------------------------------------------------------------------
  Division::Division(repr_t const & repr):
    Class("Division", repr),
    d_(arg(0), arg(1)),
    repr_(*this, d_.first, d_.second),
    str_(d_.first.str() + " / " + d_.second.str())
  {
  }

//-----------------------------------------------------------------------------
  Division::~Division()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression, Expression> const& Division::operands() const
  {
    return d_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Division::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Division::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Division::display() const
  {
  }

//-----------------------------------------------------------------------------
  Power::Power(Expression const& a, Expression const& b) :
    Class("Power"),
    apowb_(a, b),
    repr_(*this, apowb_.first, apowb_.second),
    str_(apowb_.first.str() + " ** " + apowb_.second.str())
  {
  }

//-----------------------------------------------------------------------------
  Power::Power(repr_t const & repr):
    Class("Power", repr),
    apowb_(arg(0), arg(1)),
    repr_(*this, apowb_.first, apowb_.second),
    str_(apowb_.first.str() + " ** " + apowb_.second.str())
  {
  }

//-----------------------------------------------------------------------------
  Power::~Power()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression, Expression> const& Power::operands() const
  {
    return apowb_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Power::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Power::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Power::display() const
  {
  }

//-----------------------------------------------------------------------------
  Abs::Abs(Expression const& a) :
    Class("Abs"),
    a_(a),
    repr_(*this, a_),
    str_("| " + a_.str() + " |")
  {
  }

//-----------------------------------------------------------------------------
  Abs::Abs(repr_t const & repr):
    Class("Abs", repr),
    a_(arg(0)),
    repr_(*this, a_),
    str_("| " + a_.str() + " |")
  {
  }

//-----------------------------------------------------------------------------
  Abs::~Abs()
  {
  }
  
//-----------------------------------------------------------------------------
  Expression const& Abs::operands() const
  {
    return a_;  
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Abs::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Abs::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Abs::display() const
  {
  }
}
