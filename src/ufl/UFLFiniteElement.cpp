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
                             dolfin::uint const degree) :
    FiniteElementBase("FiniteElement"),
    family_(family),
    cell_(cell),
    degree_(degree),
    value_shape_(family_.value_rank(), family_.value_rank() * cell_.geometric_dimension()),
    symmetry_(),
    sub_elements_(),
    repr_(*this, this->family(), cell, this->degree(),
          this->quadrature_scheme())
{
  // Check finite element definition
  if (!this->family().has_valid_definition(cell.domain().type(), degree))
  {
    error("The finite element definition is not valid.");
  }

  std::stringstream ssstr;
  ssstr << "<" << this->family().short_name() << degree
      << quadrature_scheme().str() << " on a " << cell.str() << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
FiniteElement::FiniteElement(repr_t const& repr) :
    FiniteElementBase("FiniteElement", repr),
    family_(arg(0)),
    cell_(arg(1)),
    degree_(arg(2)),
    value_shape_(family_.value_rank(), family_.value_rank() * cell_.geometric_dimension()),
    symmetry_(),
    sub_elements_(),
    repr_(repr)
{
  // Check finite element definition
  if (!this->family().has_valid_definition(cell().domain().type(), degree()))
  {
    error("The finite element definition is not valid.");
  }

  std::stringstream ssstr;
  ssstr << "<" << this->family().short_name() << degree()
      << quadrature_scheme().str() << " on a " << cell().str() << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
FiniteElement::~FiniteElement()
{
}

//-----------------------------------------------------------------------------
Family const& FiniteElement::family() const
{
  return family_;
}

//-----------------------------------------------------------------------------
Family::Type FiniteElement::metatype() const
{
  return family_.type();
}

//-----------------------------------------------------------------------------
Cell const& FiniteElement::cell() const
{
  return cell_;
}

//-----------------------------------------------------------------------------
type<dolfin::uint> const& FiniteElement::degree() const
{
  return degree_;
}

//-----------------------------------------------------------------------------
ValueArray const& FiniteElement::value_shape() const
{
  return value_shape_;
}

//-----------------------------------------------------------------------------
bool const FiniteElement::is_cellwise_constant() const
{
  return (family().type() == Family::R && degree() == 0);
}

//-----------------------------------------------------------------------------
std::map<dolfin::uint, dolfin::uint> const FiniteElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const FiniteElement::extract_subelement_component(
    ValueArray const& i) const
{
  return std::pair<ValueArray, ValueArray>();
}

//-----------------------------------------------------------------------------
std::pair<dolfin::uint, FiniteElementBase const *> const FiniteElement::extract_component(
    ValueArray const& i) const
{
  check_component(i);
  return std::pair<dolfin::uint, FiniteElementBase const *>(i[0], this);
}

//-----------------------------------------------------------------------------
dolfin::uint const FiniteElement::num_sub_elements() const
{
  return 0;
}

//-----------------------------------------------------------------------------
FiniteElementBase::List const& FiniteElement::sub_elements() const
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

