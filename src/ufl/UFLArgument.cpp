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
    count_(c)
  {
    std::stringstream ssrepr;
    ssrepr << "Argument("<< finite_element_.repr() << ", " << count_ << ")";
    repr_ = ssrepr.str();

    std::stringstream ssstr;
    if(count_ < 10)
    {
      ssstr  << "v_" << count_;
    }
    else
    {
      ssstr  << "v_{" << count_ << "}";
    }
    str_ = ssstr.str();
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
  bool const Argument::is_cellwise_constant() const
  {
    return true;
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

//-----------------------------------------------------------------------------
  Argument const* Argument::create(repr_t const & repr) const
  {
    if(repr.length() == 0)
      dolfin_assert("An empty signature was passed to create an Argument.");

    FiniteElementBase const * fe;
    dolfin::uint count = 0;

    std::string::const_iterator it;
    std::string help_string;
    dolfin::uint i = 0;
    for(it = repr.begin(); it!=repr.end(); ++it, ++i)
    {
      help_string += *it;
      if(help_string == "Argument(")
      {
        help_string.clear();
        std::string::const_iterator jt;
        dolfin::uint pos_end_summand = 0;
        dolfin::uint open_parentheses = 0;
        dolfin::uint close_parentheses = 0;

        for(jt = it; jt!=repr.end(); ++jt, ++pos_end_summand)
        {
          if(*jt == '(')
            open_parentheses++;
        
          if(*jt == ')')
            close_parentheses++;

          if(close_parentheses>0 && open_parentheses == close_parentheses + 1)
            break;
          std::cout << "Anzahl offene Klammern = " << open_parentheses << std::endl;
          std::cout << "Anzahl geschlossene Klammern = " << close_parentheses << std::endl;
          std::cout << "  " << pos_end_summand;
        }
        std::cout << std::endl;

        std::string string_fe = repr.substr(i+1, pos_end_summand);
        std::string string_count = repr.substr(i+pos_end_summand+3, repr.length()-pos_end_summand-i-4);
        std::cout << "substring FE  " << string_fe << std::endl;
        std::cout << "substring INDEX  " << string_count << std::endl;

        std::cout << "create FE" << std::endl;
//        fe = fe->create(string_fe);
        std::cout << "create FE done" << std::endl;
        std::cout << "create Index" << std::endl;
        count = 0;
        std::cout << "create Index done" << std::endl;
      }
    }
    return new Argument(*fe, count);
  }
}
