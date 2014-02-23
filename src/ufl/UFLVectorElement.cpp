// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLVectorElement.h>

namespace ufl
{

using dolfin::error;

//-----------------------------------------------------------------------------
VectorElement::VectorElement(Family::Type family, Cell const& cell,
                             dolfin::uint const degree, dolfin::uint const dim) :
    FiniteElementBase("VectorElement", Family::Vector, cell, degree),
    sub_element_(family, cell, degree),
    dim_(dim),
    value_shape_(),
    symmetry_(),
    sub_elements_(dim, &sub_element_),
    repr_(*this, sub_element_.family(), cell, sub_element_.degree(), dim_,
          sub_element_.quadrature_scheme())
{
  //TODO: Check mixed finite element definition

  std::stringstream ssstr;
  ssstr << "<" << sub_element_.family().short_name()
      << " vector element of degree " << degree << " on a " << cell.str()
      << ": " << sub_elements_.size() << " x " << sub_element_.str() << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
VectorElement::VectorElement(repr_t const& repr) :
    FiniteElementBase("VectorElement", repr),
    sub_element_(Family(arg(0)).type(), Cell(arg(1)), type<dolfin::uint>(arg(2))),
    dim_(arg(3)),
    value_shape_(),
    symmetry_(),
    sub_elements_(dim_, &sub_element_),
    repr_(repr)
{

  std::stringstream ssstr;
  ssstr << "<" << sub_element_.family().short_name()
      << " vector element of degree " << degree() << " on a " << cell().str()
      << ": " << sub_elements_.size() << " x " << sub_element_.str() << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
VectorElement::~VectorElement()
{
}

//-----------------------------------------------------------------------------
bool const VectorElement::is_cellwise_constant() const
{
  bool ret = true;
  for (FiniteElementBaseList::const_iterator it = sub_elements_.begin();
      it != sub_elements_.end(); ++it)
  {
    ret |= (*it)->is_cellwise_constant();
  }
  return ret;
}

//-----------------------------------------------------------------------------
std::map<dolfin::uint, dolfin::uint> const VectorElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const VectorElement::extract_subelement_component(
    ValueArray const& i) const
{
  return std::pair<dolfin::uint, dolfin::uint>();
}

//-----------------------------------------------------------------------------
std::pair<dolfin::uint, FiniteElementBase const * const> const VectorElement::extract_component(
    ValueArray const& i) const
{
  return sub_element_.extract_component(i);
}

//-----------------------------------------------------------------------------
dolfin::uint const VectorElement::num_sub_elements() const
{
  return 0;
}

//-----------------------------------------------------------------------------
FiniteElementBase::FiniteElementBaseList const& VectorElement::sub_elements() const
{
  return sub_elements_;
}

//-----------------------------------------------------------------------------
Object::repr_t const VectorElement::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const VectorElement::str() const
{
  return str_;
}

}

