// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLTensorElement.h>

namespace ufl
{

//-----------------------------------------------------------------------------
TensorElement::TensorElement(Family::Type family,
                                   Cell const& cell, uint const degree,
                                   uint const dim) :
    FiniteElementBase("TensorElement",Family::Tensor, cell, degree),
    sub_element_(family, cell, degree),
    sub_elements_(dim, &sub_element_)
{
  // Check mixed finite element definition

  // Create string representation
  QuadratureScheme qs;

  std::stringstream ssrepr;
  ssrepr << "TensorElement(" << this->family().repr()
         << ", " << cell.repr() << ", " << degree << ", " << qs.repr() << ")";
  repr_ = ssrepr.str();

  std::stringstream ssstr;
  ssstr << "<" << this->family().short_name()
        << " vector element of degree " << degree << " on a " << cell.str()
        << ": " << sub_elements_.size() << " x " << sub_element_.str() << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
TensorElement::~TensorElement()
{
}

//-----------------------------------------------------------------------------
bool const TensorElement::is_cellwise_constant() const
{
  bool ret = true;
  for ( FiniteElementBaseList::const_iterator it = sub_elements_.begin();
        it != sub_elements_.end(); ++it )
  {
    ret |= (*it)->is_cellwise_constant();
  }
  return ret;
}

//-----------------------------------------------------------------------------
std::map<dolfin::uint, dolfin::uint> const TensorElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const TensorElement::extract_subelement_component(
    ValueArray const& i) const
{
  return std::pair<dolfin::uint, dolfin::uint>();
}

//-----------------------------------------------------------------------------
std::pair<dolfin::uint, FiniteElementBase const * const> const TensorElement::extract_component(ValueArray const& i) const
{
  return sub_element_.extract_component(i);
}

//-----------------------------------------------------------------------------
dolfin::uint const TensorElement::num_sub_elements() const
{
  return 0;
}

//-----------------------------------------------------------------------------
FiniteElementBase::FiniteElementBaseList const& TensorElement::sub_elements() const
{
  return sub_elements_;
}

//-----------------------------------------------------------------------------
Object::repr_t const TensorElement::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const TensorElement::str() const
{
  return str_;
}

}

