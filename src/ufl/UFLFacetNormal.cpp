/*
 * UFLFacetNormal.cpp
 *
 *  Created on: Jan 24, 2014
 *      Author: larcher
 */

#include <dolfin/ufl/UFLFacetNormal.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLFacetNormal::UFLFacetNormal(UFLCell const& cell) :
    UFLGeometricQuantity(cell),
    shape_((cell.geometric_dimension() == 1 ? 0 : 1),cell.geometric_dimension()),
    repr_("FacetNormal(" + cell.repr() + ")"),
    str_("n")

{
}

//-----------------------------------------------------------------------------
UFLFacetNormal::~UFLFacetNormal()
{
}

//-----------------------------------------------------------------------------
std::vector<uint> const& UFLFacetNormal::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
std::string const UFLFacetNormal::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const UFLFacetNormal::str() const
{
  return str_;
}

} /* namespace dolfin */
