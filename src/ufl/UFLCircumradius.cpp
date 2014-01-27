// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLCircumradius.h>

#include <dolfin/ufl/UFLCell.h>

namespace ufl
{

//-----------------------------------------------------------------------------
Circumradius::Circumradius(Cell const& cell) :
    GeometricQuantity("Circumradius", cell),
    shape_(),
    repr_("Circumradius(" + cell.repr() + ")"),
    str_("circumradius")

{
}

//-----------------------------------------------------------------------------
Circumradius::~Circumradius()
{
}

//-----------------------------------------------------------------------------
ValueArray const& Circumradius::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
Object::repr_t const Circumradius::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const Circumradius::str() const
{
  return str_;
}

} /* namespace ufl */
