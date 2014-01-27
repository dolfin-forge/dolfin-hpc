// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLRestrictedElement.h>

namespace ufl
{

//-----------------------------------------------------------------------------
RestrictedElement::RestrictedElement(FiniteElementBase const& element,
                                           Domain::Type const domain) :
    FiniteElementBase("RestrictedElement",ElementList::Restricted,
                         element.cell(),
                         element.degree(), element.quadrature_scheme(),
                         element.value_shape()),
    element_(element)
{
  // Check mixed finite element definition

  std::stringstream ssrepr;
  ssrepr << "RestrictedElement(" << element_.repr() << ", "
         << Domain::str(domain) << ")";
  repr_ = ssrepr.str();

  std::stringstream ssstr;
  ssstr << "<" << element_.str() << ">|_" << Domain::str(domain) << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
RestrictedElement::~RestrictedElement()
{
}

//-----------------------------------------------------------------------------
bool const RestrictedElement::is_cellwise_constant() const
{
  return element_.is_cellwise_constant();
}

//-----------------------------------------------------------------------------
std::map<uint, uint> const RestrictedElement::symmetry() const
{
  return element_.symmetry();
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const RestrictedElement::extract_subelement_component(
    ValueArray const& i) const
{
  return element_.extract_subelement_component(i);
}

//-----------------------------------------------------------------------------
std::pair<uint, FiniteElementBase const * const> const RestrictedElement::extract_component(ValueArray const& i) const
{
  return element_.extract_component(i);
}

//-----------------------------------------------------------------------------
uint const RestrictedElement::num_sub_elements() const
{
  return element_.num_sub_elements();
}

//-----------------------------------------------------------------------------
FiniteElementBase::FiniteElementBaseList const& RestrictedElement::sub_elements() const
{
  return element_.sub_elements();
}

//-----------------------------------------------------------------------------
FiniteElementBase const& RestrictedElement::element()
{
  return element_;
}

//-----------------------------------------------------------------------------
Object::repr_t const RestrictedElement::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const RestrictedElement::str() const
{
  return str_;
}

}

