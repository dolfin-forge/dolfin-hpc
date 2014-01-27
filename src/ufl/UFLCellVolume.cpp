/*
 * UFLCellVolume.cpp
 *
 *  Created on: Jan 24, 2014
 *      Author: larcher
 */

#include <dolfin/ufl/UFLCellVolume.h>

#include <dolfin/ufl/UFLCell.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLCellVolume::UFLCellVolume(UFLCell const& cell) :
    UFLGeometricQuantity(cell),
    shape_(),
    repr_("CellVolume(" + cell.repr() + ")"),
    str_("volume")

{
}

//-----------------------------------------------------------------------------
UFLCellVolume::~UFLCellVolume()
{
}

//-----------------------------------------------------------------------------
ValueArray const& UFLCellVolume::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
std::string const UFLCellVolume::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const UFLCellVolume::str() const
{
  return str_;
}

} /* namespace dolfin */
