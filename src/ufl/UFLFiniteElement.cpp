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
    UFLFiniteElementBase(family, cell, degree),
    value_shape_(),
    symmetry_(),
    sub_elements_()
{
  // Check finite element definition
  if(UFLElementList::Supported().has_valid_definition(family,
                                                      cell.domain(), degree))
  {
    error("The finite element definition is not valid.");
  };

  UFLQuadratureScheme qs = "None";

  std::stringstream ssrepr;
  ssrepr << "FiniteElement("<< UFLElementList::Supported().repr(family)
         << ", "<< cell.repr() << ", " << degree << ", " << qs << ")";
  repr_ = ssrepr.str();

  std::stringstream ssstr;
  ssstr << "<" << UFLElementList::Supported().short_name(family) << degree
        << qs << " on a " << cell.repr() << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
UFLFiniteElement::~UFLFiniteElement()
{
}

//-----------------------------------------------------------------------------
bool const UFLFiniteElement::is_cellwise_constant() const
{
  return ( family() == UFLElementList::R && degree() == 0 );
}

//-----------------------------------------------------------------------------
std::map<uint, uint> const UFLFiniteElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<uint, uint> const UFLFiniteElement::extract_subelement_component(
      uint i) const
{
  return std::pair<uint, uint>();
}

//-----------------------------------------------------------------------------
uint const UFLFiniteElement::extract_component(uint i) const
{
  return i;
}

//-----------------------------------------------------------------------------
uint const UFLFiniteElement::num_sub_elements() const
{
  return 0;
}

//-----------------------------------------------------------------------------
std::vector<UFLFiniteElementBase const *> const& UFLFiniteElement::sub_elements() const
{
  return sub_elements_;
}

}

