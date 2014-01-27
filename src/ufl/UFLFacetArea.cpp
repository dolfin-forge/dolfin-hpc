/*
 * UFLFacetArea.cpp
 *
 *  Created on: Jan 24, 2014
 *      Author: larcher
 */

#include <dolfin/ufl/UFLFacetArea.h>

#include <dolfin/ufl/UFLCell.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLFacetArea::UFLFacetArea(UFLCell const& cell) :
    UFLGeometricQuantity(cell),
    shape_(),
    repr_("FacetArea(" + cell.repr() + ")"),
    str_("facetarea")

{
}

//-----------------------------------------------------------------------------
UFLFacetArea::~UFLFacetArea()
{
}

//-----------------------------------------------------------------------------
ValueArray const& UFLFacetArea::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
std::string const UFLFacetArea::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const UFLFacetArea::str() const
{
  return str_;
}

} /* namespace dolfin */
