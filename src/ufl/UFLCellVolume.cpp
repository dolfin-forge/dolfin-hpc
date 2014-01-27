// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLCellVolume.h>

#include <dolfin/ufl/UFLCell.h>

namespace ufl
{

//-----------------------------------------------------------------------------
CellVolume::CellVolume(Cell const& cell) :
    GeometricQuantity("CellVolume", cell),
    shape_(),
    repr_("CellVolume(" + cell.repr() + ")"),
    str_("volume")

{
}

//-----------------------------------------------------------------------------
CellVolume::~CellVolume()
{
}

//-----------------------------------------------------------------------------
ValueArray const& CellVolume::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
Object::repr_t const CellVolume::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const CellVolume::str() const
{
  return str_;
}

} /* namespace ufl */
