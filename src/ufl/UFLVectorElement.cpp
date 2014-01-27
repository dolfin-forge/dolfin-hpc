// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLVectorElement.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLVectorElement::UFLVectorElement(UFLElementList::FamilyType family,
                                   UFLCell const& cell, uint const degree,
                                   uint const dim) :
    UFLFiniteElementBase(UFLElementList::Vector, cell, degree),
    sub_element_(family, cell, degree),
    sub_elements_(dim, &sub_element_)
{
  // Check mixed finite element definition

  // Create string representation
  UFLQuadratureScheme qs = "None";

  std::stringstream ssrepr;
  ssrepr << "VectorElement(" << UFLElementList::Supported().repr(family)
         << ", " << cell.repr() << ", " << degree << ", " << qs << ")";
  repr_ = ssrepr.str();

  std::stringstream ssstr;
  ssstr << "<" << UFLElementList::Supported().short_name(family)
        << " vector element of degree " << degree << " on a " << cell.str()
        << ": " << sub_elements_.size() << " x " << sub_element_.str() << ">";
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
UFLVectorElement::~UFLVectorElement()
{
}

//-----------------------------------------------------------------------------
bool const UFLVectorElement::is_cellwise_constant() const
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
std::map<uint, uint> const UFLVectorElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const UFLVectorElement::extract_subelement_component(
    ValueArray const& i) const
{
  return std::pair<uint, uint>();
}

//-----------------------------------------------------------------------------
std::pair<uint, UFLFiniteElementBase const * const> const UFLVectorElement::extract_component(ValueArray const& i) const
{
  return sub_element_.extract_component(i);
}

//-----------------------------------------------------------------------------
uint const UFLVectorElement::num_sub_elements() const
{
  return 0;
}

//-----------------------------------------------------------------------------
UFLFiniteElementBase::FiniteElementBaseList const& UFLVectorElement::sub_elements() const
{
  return sub_elements_;
}

}

