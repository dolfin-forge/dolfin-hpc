/*
 * UFLSpatialCoordinate.cpp
 *
 *  Created on: Jan 24, 2014
 *      Author: larcher
 */

#include <dolfin/ufl/UFLSpatialCoordinate.h>

#include <dolfin/ufl/UFLCell.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLSpatialCoordinate::UFLSpatialCoordinate(UFLCell const& cell) :
    UFLGeometricQuantity(cell),
    shape_((cell.geometric_dimension() == 1 ? 0 : 1),cell.geometric_dimension()),
    repr_("SpatialCoordinate(" + cell.repr() + ")"),
    str_("x")
{
}

//-----------------------------------------------------------------------------
UFLSpatialCoordinate::~UFLSpatialCoordinate()
{
}

//-----------------------------------------------------------------------------
std::vector<uint> const& UFLSpatialCoordinate::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
std::string const UFLSpatialCoordinate::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const UFLSpatialCoordinate::str() const
{
  return str_;
}

} /* namespace dolfin */
