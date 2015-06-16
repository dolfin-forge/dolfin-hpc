// Copyright (C) 2014 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#include <dolfin/ufl/UFLList.h>

namespace ufl
{

//-----------------------------------------------------------------------------
  List::List(std::vector<Integral const *> const& integrals) :
    Class("[", "]"),
    integrals_(integrals),
    repr_(*this, integrals_),
    str_("")
  {
  }

//-----------------------------------------------------------------------------
  List::List(repr_t const& repr) :
    Class("[", "]", repr),
    integrals_(fill_expressions(args())),
    repr_(*this, integrals_),
    str_("")
  {
  }

//-----------------------------------------------------------------------------
  List::~List()
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Integral const *> const& List::get_integrals() const
  {
    return integrals_;
  }

//-----------------------------------------------------------------------------
  Object::repr_t const& List::repr() const
  {
    return repr_;
  }

  //-----------------------------------------------------------------------------
  std::string const& List::str() const
  {
    return str_;
  }

//-----------------------------------------------------------------------------
  void List::display() const
  {
  }

//-----------------------------------------------------------------------------
  std::vector<Integral const *> const List::fill_expressions(std::vector<repr_t> const& reprs)
  {
    std::vector<Integral const *> expr;
    for(dolfin::uint i=0; i<reprs.size(); ++i)
      expr.push_back(new Integral(reprs[i]));
    return expr;
  }
}
