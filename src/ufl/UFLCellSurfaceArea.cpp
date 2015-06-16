// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLCellSurfaceArea.h>

#include <dolfin/ufl/UFLCell.h>

namespace ufl
{

//-----------------------------------------------------------------------------
CellSurfaceArea::CellSurfaceArea(Cell const& cell) :
    GeometricQuantity("CellSurfaceArea", cell),
    shape_(),
    repr_(*this, cell),
    str_("surfacearea")

{
}

//-----------------------------------------------------------------------------
CellSurfaceArea::~CellSurfaceArea()
{
}

//-----------------------------------------------------------------------------
ValueArray const& CellSurfaceArea::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
std::string const& CellSurfaceArea::str() const
{
  return str_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& CellSurfaceArea::repr() const
{
  return repr_;
}

} /* namespace ufl */
