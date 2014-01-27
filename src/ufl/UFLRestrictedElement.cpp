// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLRestrictedElement.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLRestrictedElement::UFLRestrictedElement(UFLFiniteElementBase const& element,
                                           UFLDomain::Type const domain) :
    UFLFiniteElementBase(UFLElementList::Restricted, element.cell(),
                         element.degree(), element.quadrature_scheme(),
                         element.value_shape()),
    element_(element)
{
  // Check mixed finite element definition

  UFLQuadratureScheme qs = "None";

  std::stringstream ssrepr;
  ssrepr << "RestrictedElement(" << element_.repr() << ", "
         << UFLDomain::str(domain) << ")";
  repr_ = ssrepr.str();

  std::stringstream ssstr;
  ssstr << "<" << element_.str() << ">|_" << UFLDomain::str(domain) << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
UFLRestrictedElement::~UFLRestrictedElement()
{
}

//-----------------------------------------------------------------------------
bool const UFLRestrictedElement::is_cellwise_constant() const
{
  return element_.is_cellwise_constant();
}

//-----------------------------------------------------------------------------
std::map<uint, uint> const UFLRestrictedElement::symmetry() const
{
  return element_.symmetry();
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const UFLRestrictedElement::extract_subelement_component(
    ValueArray const& i) const
{
  return element_.extract_subelement_component(i);
}

//-----------------------------------------------------------------------------
std::pair<uint, UFLFiniteElementBase const * const> const UFLRestrictedElement::extract_component(ValueArray const& i) const
{
  return element_.extract_component(i);
}

//-----------------------------------------------------------------------------
uint const UFLRestrictedElement::num_sub_elements() const
{
  return element_.num_sub_elements();
}

//-----------------------------------------------------------------------------
UFLFiniteElementBase::FiniteElementBaseList const& UFLRestrictedElement::sub_elements() const
{
  return element_.sub_elements();
}

//-----------------------------------------------------------------------------
UFLFiniteElementBase const& UFLRestrictedElement::element()
{
  return element_;
}

}

