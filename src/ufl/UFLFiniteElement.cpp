    initeElementBase("FiniteElement", family, cell, degree),
// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLFiniteElement.h>

namespace ufl
{

using dolfin::error;

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(Family::Type family, Cell const& cell,
                             uint const degree) :
    FiniteElementBase("FiniteElement", family, cell, degree),
    value_shape_(),
    symmetry_(),
    sub_elements_()
{
  // Check finite element definition
  if (!this->family().has_valid_definition(cell.domain().type(), degree))
  {
    error("The finite element definition is not valid.");
  }

  QuadratureScheme qs;

  std::stringstream ssrepr;
  ssrepr << "FiniteElement(" << this->family().repr() << ", " << cell.repr()
      << ", " << degree << ", " << qs.repr() << ")";
  repr_ = ssrepr.str();

  std::stringstream ssstr;
  ssstr << "<" << this->family().short_name() << degree << qs.str() << " on a "
      << cell.str() << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(repr_t const& repr) :
    FiniteElementBase("FiniteElement", repr),
    value_shape_(),
    symmetry_(),
    sub_elements_()
{
  // Check finite element definition
  if (!this->family().has_valid_definition(cell().domain().type(), degree()))
  {
    error("The finite element definition is not valid.");
  }

  QuadratureScheme qs;

  repr_ = repr;

  std::stringstream ssstr;
  ssstr << "<" << this->family().short_name() << degree() << qs.str() << " on a "
      << cell().str() << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
FiniteElement::~FiniteElement()
{
}

//-----------------------------------------------------------------------------
bool const FiniteElement::is_cellwise_constant() const
{
  return (family().type() == Family::R && degree() == 0);
}

//-----------------------------------------------------------------------------
std::map<uint, uint> const FiniteElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const FiniteElement::extract_subelement_component(
    ValueArray const& i) const
{
  return std::pair<uint, uint>();
}

//-----------------------------------------------------------------------------
std::pair<uint, FiniteElementBase const * const > const FiniteElement::extract_component(
    ValueArray const& i) const
{
  check_component(i);
  return std::pair<uint, FiniteElementBase const * const >(i[0], this);
}

//-----------------------------------------------------------------------------
uint const FiniteElement::num_sub_elements() const
{
  return 0;
}

//-----------------------------------------------------------------------------
FiniteElementBase::FiniteElementBaseList const& FiniteElement::sub_elements() const
{
  return sub_elements_;
}

//-----------------------------------------------------------------------------
Object::repr_t const FiniteElement::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const FiniteElement::str() const
{
  return str_;
}

}

