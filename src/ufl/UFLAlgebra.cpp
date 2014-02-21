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
    s1_(s1),
    s2_(s2),
    repr_(*this, s1_, s2_),
    str_(s1_.str() + " + " + s2_.str())
  {
  }

//-----------------------------------------------------------------------------
  Sum::Sum(repr_t const & repr):
    Class("Sum", repr),
    s1_(arg(0)),
    s2_(arg(1)),
    repr_(*this, s1_, s2_),
    str_(s1_.str() + " + " + s2_.str())
  {
  }

//-----------------------------------------------------------------------------
  Sum::~Sum()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression const, Expression const> const& Sum::operands() const
  {
    return std::make_pair(s1_, s2_);  
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
    p1_(p1),
    p2_(p2),
    repr_(*this, p1_, p2_),
    str_(p1_.str() + " * " + p2_.str())
  {
  }

//-----------------------------------------------------------------------------
  Product::Product(repr_t const & repr):
    Class("Product", repr),
    p1_(arg(0)),
    p2_(arg(1)),
    repr_(*this, p1_, p2_),
    str_(p1_.str() + " * " + p2_.str())
  {
  }

//-----------------------------------------------------------------------------
  Product::~Product()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression const, Expression const> const& Product::operands() const
  {
    return std::make_pair(p1_, p2_);  
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
    d1_(d1),
    d2_(d2),
    repr_(*this, d1_, d2_),
    str_(d1_.str() + " / " + d2_.str())
  {
  }

//-----------------------------------------------------------------------------
  Division::Division(repr_t const & repr):
    Class("Division", repr),
    d1_(arg(0)),
    d2_(arg(1)),
    repr_(*this, d1_, d2_),
    str_(d1_.str() + " / " + d2_.str())
  {
  }

//-----------------------------------------------------------------------------
  Division::~Division()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression const, Expression const> const& Division::operands() const
  {
    return std::make_pair(d1_, d2_);  
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
    a_(a),
    b_(b),
    repr_(*this, a_, b_),
    str_(a_.str() + " ** " + b_.str())
  {
  }

//-----------------------------------------------------------------------------
  Power::Power(repr_t const & repr):
    Class("Power", repr),
    a_(arg(0)),
    b_(arg(1)),
    repr_(*this, a_, b_),
    str_(a_.str() + " ** " + b_.str())
  {
  }

//-----------------------------------------------------------------------------
  Power::~Power()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression const, Expression const> const& Power::operands() const
  {
    return std::make_pair(a_, b_);  
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
