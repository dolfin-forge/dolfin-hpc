// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLIndexSum.h>

//#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  IndexSum::IndexSum(Expression const& summand, IndexBase const& index) :
    Class("Form"),
    summand_(summand),
    index_(index)
  {
    std::stringstream ssrepr;
    ssrepr << "IndexSum(" << summand_.repr() << "," << index_.repr() << ")";
    repr_ = ssrepr.str();

    //Is this the same implementation as in python?
    std::stringstream ssstr;
    ssstr << "sum_{" << index_.str() << "} " << summand_.str() << " ";
    str_ = ssstr.str();
  }

//-----------------------------------------------------------------------------
  IndexSum::~IndexSum()
  {
  }
  
//-----------------------------------------------------------------------------
  IndexBase const& IndexSum::index() const
  {
    return index_;  
  }

//-----------------------------------------------------------------------------
  dolfin::uint const& IndexSum::dimension() const
  {
    return 0;  
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression const, IndexBase const> const& IndexSum::operands() const
  {
    return std::make_pair(summand_, index_);  
  }

//-----------------------------------------------------------------------------
  Object::repr_t const IndexSum::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const IndexSum::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void IndexSum::display() const
  {
  }


//-----------------------------------------------------------------------------
  IndexSum const* IndexSum::create(repr_t const & repr) const
  {
    if(repr.length() == 0)
      dolfin_assert("An empty signature was passed to create an IndexSum.");

    Expression const * summand;
    Index const * index;

    std::string::const_iterator it;
    std::string help_string;
    dolfin::uint i = 0;
    for(it = repr.begin(); it!=repr.end(); ++it, ++i)
    {
      help_string += *it;
      if(help_string == "IndexSum(")
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
//          std::cout << "Anzahl offene Klammern = " << open_parentheses << std::endl;
//          std::cout << "Anzahl geschlossene Klammern = " << close_parentheses << std::endl;
//          std::cout << "  " << pos_end_summand;
        }
        std::cout << std::endl;

        std::string string_summand = repr.substr(i+1, pos_end_summand);
        std::string string_index = repr.substr(i+pos_end_summand+3, repr.length()-pos_end_summand-i-4);
        std::cout << "substring SUMMAND  " << string_summand << std::endl;
        std::cout << "substring INDEX  " << string_index << std::endl;

        std::cout << "create Expression" << std::endl;
        summand = summand->create(string_summand);
        std::cout << "create Expression done" << std::endl;
        std::cout << "create Index" << std::endl;
        index = index->create(string_index);
        std::cout << "create Index done" << std::endl;
      }
    }
    return new IndexSum(*summand, *index);
  }
}
