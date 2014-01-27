// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLFacetNormal.h>

namespace ufl
{

//-----------------------------------------------------------------------------
FacetNormal::FacetNormal(Cell const& cell) :
    GeometricQuantity(cell),
    shape_((cell.geometric_dimension() == 1 ? 0 : 1),cell.geometric_dimension()),
    repr_("FacetNormal(" + cell.repr() + ")"),
    str_("n")

{
}

//-----------------------------------------------------------------------------
FacetNormal::~FacetNormal()
{
}

//-----------------------------------------------------------------------------
ValueArray const& FacetNormal::shape() const
{
  return shape_;
}

//-----------------------------------------------------------------------------
std::string const FacetNormal::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const FacetNormal::str() const
{
  return str_;
}

} /* namespace ufl */
