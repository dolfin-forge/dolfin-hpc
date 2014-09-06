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
    summand_mindex_(summand, index)
  {
    std::stringstream ssrepr;
    ssrepr << "IndexSum(" << summand.repr() << "," << index.repr() << ")";
    repr_ = ssrepr.str();

    //Is this the same implementation as in python?
    std::stringstream ssstr;
    ssstr << "sum_{" << index.str() << "} " << summand.str() << " ";
    str_ = ssstr.str();
  }

//-----------------------------------------------------------------------------
  IndexSum::IndexSum(repr_t const & repr) :
    Class("IndexSum", repr),
    summand_mindex_(arg(0), arg(1))
  {
  }

//-----------------------------------------------------------------------------
  IndexSum::~IndexSum()
  {
  }
  
//-----------------------------------------------------------------------------
  MultiIndex const& IndexSum::index() const
  {
    return summand_mindex_.second;
  }

//-----------------------------------------------------------------------------
  dolfin::uint IndexSum::dimension() const
  {
    return 0;  
  }
  
//-----------------------------------------------------------------------------
  std::pair<Expression, MultiIndex> const& IndexSum::operands() const
  {
    return summand_mindex_;
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
