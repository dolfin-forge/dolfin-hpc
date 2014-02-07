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
//    if(repr.length() == 0)
//      dolfin_assert("An empty signature was passed to create a Product.");
//
//    Expression const * p1;
//    Expression const * p2;
//
//    std::string::const_iterator it;
//    std::string help_string;
//    dolfin::uint i= 0;
//    for(it = repr.begin(); it!=repr.end(); ++it, ++i)
//    {
//      help_string += *it;
//      if(help_string == "Product(")
//      {
//        help_string.clear();
//        std::string::const_iterator jt;
//        dolfin::uint pos_end_summand = 0;
//        dolfin::uint open_parentheses = 0;
//        dolfin::uint close_parentheses = 0;
//
//        for(jt = it; jt!=repr.end(); ++jt, ++pos_end_summand)
//        {
//          if(*jt == '(')
//            open_parentheses++;
//        
//          if(*jt == ')')
//            close_parentheses++;
//
//          if(close_parentheses>0 && open_parentheses == close_parentheses + 1)
//            break;
//          std::cout << "Anzahl offene Klammern = " << open_parentheses << std::endl;
//          std::cout << "Anzahl geschlossene Klammern = " << close_parentheses << std::endl;
//          std::cout << "                  " << pos_end_summand << std::endl;
//        }
//
//        std::string string_p1 = repr.substr(i+1, pos_end_summand);
//        std::string string_p2 = repr.substr(i+pos_end_summand+3, repr.length()-pos_end_summand-i-4);
//        std::cout << "substring P1  " << string_p1 << std::endl;
//        std::cout << "substring P2  " << string_p2 << std::endl;
//        std::cout << "create Expression 1" << std::endl;
//        p1 = p1->create(string_p1);
//        std::cout << "create Expression 1 done" << std::endl;
//        std::cout << "create Expression 2" << std::endl;
//        p2 = p2->create(string_p2);
//        std::cout << "create Expression 2 done" << std::endl;
//      }
//    }
//    return new Product(*p1, *p2);
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


}
