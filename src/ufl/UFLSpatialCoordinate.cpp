// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLSpatialCoordinate.h>

#include <dolfin/ufl/UFLCell.h>

namespace ufl
{

//-----------------------------------------------------------------------------
SpatialCoordinate::SpatialCoordinate(Cell const& cell) :
    GeometricQuantity("SpatialCoordinate", cell),
    shape_((cell.geometric_dimension() == 1 ? 0 : 1),cell.geometric_dimension()),
    repr_(*this, cell),
    str_("x")
{
}

//-----------------------------------------------------------------------------
SpatialCoordinate::~SpatialCoordinate()
{
}

//-----------------------------------------------------------------------------
ValueArray const& SpatialCoordinate::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& SpatialCoordinate::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& SpatialCoordinate::str() const
{
  return str_;
}

} /* namespace ufl */
