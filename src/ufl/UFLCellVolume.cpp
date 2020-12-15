// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLCellVolume.h>

#include <dolfin/ufl/UFLCell.h>

namespace ufl
{

//-----------------------------------------------------------------------------
CellVolume::CellVolume(Cell const& cell) :
    GeometricQuantity("CellVolume", cell),
    shape_(),
    repr_(*this, cell),
    str_("volume")

{
}

//-----------------------------------------------------------------------------
CellVolume::~CellVolume()
= default;

//-----------------------------------------------------------------------------
ValueArray const& CellVolume::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& CellVolume::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& CellVolume::str() const
{
  return str_;
}

} /* namespace ufl */
