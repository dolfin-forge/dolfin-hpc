// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLFiniteElement.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLFiniteElement::UFLFiniteElement(UFLElementList::FamilyType family,
                                           UFLCell const& cell,
                                           uint const degree) :
    UFLFiniteElementBase(family, cell, degree)
{
  // Check finite element definition
  if(UFLElementList::Supported().has_valid_definition(family,
                                                      cell.domain(), degree))
  {
    error("The finite element definition is not valid.");
  };
}

//-----------------------------------------------------------------------------
UFLFiniteElement::~UFLFiniteElement()
{
}

//-----------------------------------------------------------------------------
UFLElementList::FamilyType const UFLFiniteElement::family() const
{
  return family_;
}

//-----------------------------------------------------------------------------
UFLCell const UFLFiniteElement::cell() const
{
  return cell_;
}

//-----------------------------------------------------------------------------
uint const UFLFiniteElement::degree() const
{
  return degree_;
}

//-----------------------------------------------------------------------------
ValueShape const UFLFiniteElement::value_shape() const
{
  return value_shape_;
}

}

