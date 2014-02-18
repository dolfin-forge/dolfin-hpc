// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLExpression.h>
#include <dolfin/ufl/UFLTuple.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  Tuple::Tuple(std::vector<Expression const *>& expressions) :
    Class("Tuple"),
    expressions_(expressions),
//    repr_(*this, expressions_),
    str_("")//Tuple(*(" + expressions.str() + "))")
  {
//    expressions_.clear();
//    expressions_.resize(integrals.size());
//    for(dolfin::uint i=0; i<args.size(); ++i)
//      expressions_[i] = new Expression(*args[i]);
  }

//-----------------------------------------------------------------------------
  Tuple::Tuple(repr_t const & repr) :
    Class("Tuple", repr),
//    repr_(*this, arg(0)),
    str_("")//Tuple(*(" + arg(0) + "))")
  {
//    std::cout << "C Tuple" << std::endl;
    std::vector<repr_t> const arguments = args();
    expressions_.clear();
    expressions_.resize(arguments.size());
    for(dolfin::uint i=0; i<arguments.size(); ++i)
      expressions_[i]= new Expression(arguments[i]);
//    std::cout << "C Tuple end" << std::endl;
  }
//-----------------------------------------------------------------------------
  Tuple::~Tuple()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const Tuple::repr() const
  {
    return repr_;
  }

  //-----------------------------------------------------------------------------
  std::string const Tuple::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Tuple::display() const
  {
  }
}
