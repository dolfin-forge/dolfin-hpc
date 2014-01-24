/*
 * UFLCircumradius.cpp
 *
 *  Created on: Jan 24, 2014
 *      Author: larcher
 */

#include <dolfin/ufl/UFLCircumradius.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLCircumradius::UFLCircumradius(UFLCell const& cell) :
    UFLGeometricQuantity(cell),
    shape_(),
    repr_("Circumradius(" + cell.repr() + ")"),
    str_("circumradius")

{
}

//-----------------------------------------------------------------------------
UFLCircumradius::~UFLCircumradius()
{
}

//-----------------------------------------------------------------------------
std::vector<uint> const& UFLCircumradius::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
std::string const UFLCircumradius::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const UFLCircumradius::str() const
{
  return str_;
}

} /* namespace dolfin */
