// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLTensorElement.h>

namespace ufl
{

//-----------------------------------------------------------------------------
TensorElement::TensorElement(ElementList::FamilyType family,
                                   Cell const& cell, uint const degree,
                                   uint const dim) :
    FiniteElementBase(ElementList::Tensor, cell, degree),
    sub_element_(family, cell, degree),
    sub_elements_(dim, &sub_element_)
{
  // Check mixed finite element definition

  // Create string representation
  QuadratureScheme qs;

  std::stringstream ssrepr;
  ssrepr << "TensorElement(" << ElementList::Supported().repr(family)
         << ", " << cell.repr() << ", " << degree << ", " << qs.repr() << ")";
  repr_ = ssrepr.str();

  std::stringstream ssstr;
  ssstr << "<" << ElementList::Supported().short_name(family)
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
std::map<uint, uint> const TensorElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const TensorElement::extract_subelement_component(
    ValueArray const& i) const
{
  return std::pair<uint, uint>();
}

//-----------------------------------------------------------------------------
std::pair<uint, FiniteElementBase const * const> const TensorElement::extract_component(ValueArray const& i) const
{
  return sub_element_.extract_component(i);
}

//-----------------------------------------------------------------------------
uint const TensorElement::num_sub_elements() const
{
  return 0;
}

//-----------------------------------------------------------------------------
FiniteElementBase::FiniteElementBaseList const& TensorElement::sub_elements() const
{
  return sub_elements_;
}

}

