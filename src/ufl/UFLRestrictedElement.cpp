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
                                     Domain const& domain) :
    FiniteElementBase("RestrictedElement", element.quadrature_scheme()),
    element_(element),
    family_(Family::Mixed),
    cell_(element.cell()),
    degree_(element.degree()),
    value_shape_(element.value_shape())
{
  // Check mixed finite element definition

  std::stringstream ssrepr;
  ssrepr << "RestrictedElement(" << element_.repr() << ", " << domain.repr()
      << ")";
  repr_ = ssrepr.str();

  std::stringstream ssstr;
  ssstr << "<" << element_.str() << ">|_" << domain.str() << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
RestrictedElement::~RestrictedElement()
{
}

//-----------------------------------------------------------------------------
Family const& RestrictedElement::family() const
{
  return family_;
}

//-----------------------------------------------------------------------------
Family::Type RestrictedElement::metatype() const
{
  return Family::Restricted;
}

//-----------------------------------------------------------------------------
Cell const& RestrictedElement::cell() const
{
  return cell_;
}

//-----------------------------------------------------------------------------
type<dolfin::uint> const& RestrictedElement::degree() const
{
  return degree_;
}

//-----------------------------------------------------------------------------
ValueArray const& RestrictedElement::value_shape() const
{
  return value_shape_;
}

//-----------------------------------------------------------------------------
bool RestrictedElement::is_cellwise_constant() const
{
  return element_.is_cellwise_constant();
}

//-----------------------------------------------------------------------------
std::map<dolfin::uint, dolfin::uint> const& RestrictedElement::symmetry() const
{
  return element_.symmetry();
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> RestrictedElement::extract_subelement_component(
    ValueArray const& i) const
{
  return element_.extract_subelement_component(i);
}

//-----------------------------------------------------------------------------
std::pair<dolfin::uint, FiniteElementBase const *> RestrictedElement::extract_component(
    ValueArray const& i) const
{
  return element_.extract_component(i);
}

//-----------------------------------------------------------------------------
dolfin::uint RestrictedElement::num_sub_elements() const
{
  return element_.num_sub_elements();
}

//-----------------------------------------------------------------------------
FiniteElementBase::List const& RestrictedElement::sub_elements() const
{
  return element_.sub_elements();
}

//-----------------------------------------------------------------------------
FiniteElementBase const& RestrictedElement::element()
{
  return element_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& RestrictedElement::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& RestrictedElement::str() const
{
  return str_;
}

}

