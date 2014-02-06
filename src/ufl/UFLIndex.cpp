// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLIndex.h>

//#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  IndexBase::IndexBase(dolfin::uint const& count) :
    count_(count)
  {
  }

//-----------------------------------------------------------------------------
  IndexBase::~IndexBase()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const IndexBase::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const IndexBase::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void IndexBase::display() const
  {
  }


//-----------------------------------------------------------------------------
  IndexBase const* IndexBase::create(repr_t const & repr) const
  {
    /*
    if(repr.length() == 0)
      dolfin_assert("An empty signature was passed to create an IndexSum.");

    Expression const * summand;
    Index const * index;

    std::string::const_iterator it;
    std::string help_string;
    dolfin::uint index = 0;
    for(it = repr.begin(); it!=repr.end(); ++it, ++index)
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

          if(open_parentheses>0 && open_parentheses == close_parentheses)
            break;
          std::cout << "Anzahl offene Klammern = " << open_parentheses << std::endl;
          std::cout << "Anzahl geschlossene Klammern = " << close_parentheses << std::endl;
          std::cout << "                  " << pos_end_summand << std::endl;
        }

        std::string string_summand = repr.substr(index+1, pos_end_summand);
        std::string string_index = repr.substr(pos_end_summand+1, repr.length());
        std::cout << "substring SUMMAND  " << string_summand << std::endl;
        std::cout << "substring INDEX  " << string_index << std::endl;
        summand = summand->create(string_summand);
        index = index->create(string_index);
      }
    }
    */
    return new IndexBase(0);
  }
  
//-----------------------------------------------------------------------------
  Index::Index(dolfin::uint const& count) :
    IndexBase(count)
  {
  }

//-----------------------------------------------------------------------------
  Index::~Index()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const Index::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const Index::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void Index::display() const
  {
  }


//-----------------------------------------------------------------------------
  Index const* Index::create(repr_t const & repr) const
  {
    if(repr.length() == 0)
      dolfin_assert("An empty signature was passed to create an Index.");

    std::string::const_iterator it;
    std::string help_string;
    dolfin::uint index = 0;
    for(it = repr.begin(); it!=repr.end(); ++it, ++index)
    {
      help_string += *it;
      if(help_string == "Index(")
      {
        help_string.clear();

        std::string sub_string = repr.substr(index+1, 1);
        std::cout << "substring " << sub_string << std::endl;
      }
    }
    return new Index(index);
  }

//-----------------------------------------------------------------------------
  FixedIndex::FixedIndex(dolfin::uint const& count) :
    IndexBase(count)
  {
  }

//-----------------------------------------------------------------------------
  FixedIndex::~FixedIndex()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const FixedIndex::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const FixedIndex::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void FixedIndex::display() const
  {
  }


//-----------------------------------------------------------------------------
  FixedIndex const* FixedIndex::create(repr_t const & repr) const
  {
    if(repr.length() == 0)
      dolfin_assert("An empty signature was passed to create a FixedIndex.");

    std::string::const_iterator it;
    std::string help_string;
    dolfin::uint index = 0;
    for(it = repr.begin(); it!=repr.end(); ++it, ++index)
    {
      help_string += *it;
      if(help_string == "FixedIndex(")
      {
        help_string.clear();

        std::string sub_string = repr.substr(index+1, 1);
        std::cout << "substring " << sub_string << std::endl;
      }
    }
    return new FixedIndex(index);
  }
  
//-----------------------------------------------------------------------------
  MultiIndex::MultiIndex(dolfin::uint const& count) :
    IndexBase(count)    
  {
  }

//-----------------------------------------------------------------------------
  MultiIndex::MultiIndex(IndexBase const& index) :
    IndexBase(index)
  {
  }

//-----------------------------------------------------------------------------
  MultiIndex::~MultiIndex()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const MultiIndex::repr() const
  {
    return repr_;
  }

//-----------------------------------------------------------------------------
  std::string const MultiIndex::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void MultiIndex::display() const
  {
  }


//-----------------------------------------------------------------------------
  MultiIndex const* MultiIndex::create(repr_t const & repr) const
  {
    if(repr.length() == 0)
      dolfin_assert("An empty signature was passed to create a MultiIndex.");

    std::string::const_iterator it;
    std::string help_string;
    dolfin::uint index = 0;
    for(it = repr.begin(); it!=repr.end(); ++it, ++index)
    {
      help_string += *it;
      if(help_string == "MultiIndex(")
      {
        help_string.clear();

        std::string sub_string = repr.substr(index+1, 1);
        std::cout << "substring " << sub_string << std::endl;
      }
    }
    return new MultiIndex(index);
  }
}
