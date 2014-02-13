// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed: 

#include <dolfin/ufl/UFLList.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  List::List(std::vector<Integral> const & integrals) :
    Class("[", "]")
  {
    integrals_.clear();
    integrals_.resize(integrals.size());
    for(dolfin::uint i=0; i<integrals.size(); ++i)
      integrals_[i] = new Integral(integrals[i]);
  }

//-----------------------------------------------------------------------------
  List::List(repr_t const & repr) :
    Class("[", "]", repr)
  {
//    std::cout << "C List" << std::endl;
    std::vector<repr_t> const arguments = args();
    integrals_.clear();
    integrals_.resize(arguments.size());
    for(dolfin::uint i=0; i<arguments.size(); ++i)
      integrals_[i]= new Integral(arguments[i]);
//    std::cout << "C List end" << std::endl;
  }
//-----------------------------------------------------------------------------
  List::~List()
  {
  }
  
//-----------------------------------------------------------------------------
  Object::repr_t const List::repr() const
  {
    return repr_;
  }

  //-----------------------------------------------------------------------------
  std::string const List::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void List::display() const
  {
  }
}
