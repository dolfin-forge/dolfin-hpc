// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLFacetArea.h>

#include <dolfin/ufl/UFLCell.h>

namespace ufl
{

//-----------------------------------------------------------------------------
FacetArea::FacetArea(Cell const& cell) :
    GeometricQuantity("FacetArea", cell),
    shape_(),
    repr_(*this, cell),
    str_("facetarea")

{
}

//-----------------------------------------------------------------------------
FacetArea::~FacetArea()
{
}

//-----------------------------------------------------------------------------
ValueArray const& FacetArea::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& FacetArea::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& FacetArea::str() const
{
  return str_;
}

} /* namespace ufl */
