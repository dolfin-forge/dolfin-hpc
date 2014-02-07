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
    Class("IndexSum"),
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
  IndexSum::IndexSum(repr_t const & repr) :
    Class("IndexSum", repr),
    summand_(arg(0)),
    index_(arg(1))
  {
  }

//-----------------------------------------------------------------------------
  IndexSum::~IndexSum()
  {
  }
  
//-----------------------------------------------------------------------------
  MultiIndex const& IndexSum::index() const
  {
    return index_;  
  }

//-----------------------------------------------------------------------------
  dolfin::uint const& IndexSum::dimension() const
  {
    return 0;  
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression const, MultiIndex const> const& IndexSum::operands() const
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

}
