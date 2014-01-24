/*
 * UFLCellSurfaceArea.cpp
 *
 *  Created on: Jan 24, 2014
 *      Author: larcher
 */

#include <dolfin/ufl/UFLCellSurfaceArea.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLCellSurfaceArea::UFLCellSurfaceArea(UFLCell const& cell) :
    UFLGeometricQuantity(cell),
    shape_(),
    repr_("CellSurfaceArea(" + cell.repr() + ")"),
    str_("surfacearea")

{
}

//-----------------------------------------------------------------------------
UFLCellSurfaceArea::~UFLCellSurfaceArea()
{
}

//-----------------------------------------------------------------------------
std::vector<uint> const& UFLCellSurfaceArea::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
std::string const UFLCellSurfaceArea::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const UFLCellSurfaceArea::str() const
{
  return str_;
}

} /* namespace dolfin */
