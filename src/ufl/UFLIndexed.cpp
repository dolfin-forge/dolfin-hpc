// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLIndexed.h>

//#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  Indexed::Indexed(Expression const& expression, Index const& index) :
    Class("Indexed"),
    expr_index_(expression, index)
  {
    std::stringstream ssrepr;
    ssrepr << "Indexed(" << expression.repr() << "," << index.repr() << ")";
    repr_ = ssrepr.str();

    //Is this the same implementation as in python?
    std::stringstream ssstr;
    ssstr << expression.str() << "[" << index.str() << "]";
    str_ = ssstr.str();
  }

//-----------------------------------------------------------------------------
  Indexed::Indexed(repr_t const & repr) :
    Class("Indexed", repr),
    expr_index_(arg(0), arg(1))
  {
//    if(repr.length() == 0)
//      dolfin_assert("An empty signature was passed to create an Indexed.");
//
//    Expression const * expression;
//    IndexBase const * index;
//
//    std::string::const_iterator it;
//    std::string help_string;
//    dolfin::uint i = 0;
//    for(it = repr.begin(); it!=repr.end(); ++it, ++i)
//    {
//      help_string += *it;
//      if(help_string == "Indexed(")
//      {
//        help_string.clear();
//        std::string::const_iterator jt;
//        dolfin::uint pos_end_expr = 0;
//        dolfin::uint open_parentheses = 0;
//        dolfin::uint close_parentheses = 0;
//
//        for(jt = it; jt!=repr.end(); ++jt, ++pos_end_expr)
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
//          std::cout << "  " << pos_end_expr;
//        }
//        std::cout << std::endl;
//
//        std::string string_expr = repr.substr(i+1, pos_end_expr);
//        std::string string_index = repr.substr(i+pos_end_expr+3, repr.length()-pos_end_expr-i-4);
//        std::cout << "substring EXPR  " << string_expr << std::endl;
//        std::cout << "substring INDEX  " << string_index << std::endl;
//        expression = expression->create(string_expr);
//        index = index->create(string_index);
//      }
//    }
//    return new Indexed(*expression, *index);
  }

//-----------------------------------------------------------------------------
  Indexed::~Indexed()
  {
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression, Index> const& Indexed::operands() const
  {
    return expr_index_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const Indexed::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Indexed::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Indexed::display() const
  {
  }
}
