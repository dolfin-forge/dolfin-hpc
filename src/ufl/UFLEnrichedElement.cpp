// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLEnrichedElement.h>
#include <dolfin/ufl/UFLElementList.h>

namespace ufl
{

//-----------------------------------------------------------------------------
EnrichedElement::EnrichedElement(
    FiniteElementBaseList const& elements ) :
    FiniteElementBase(ElementList::Enriched, get_cell(elements),
                         get_degree_max(elements)),
    sub_elements_(elements)
{
  // Create string representation
  std::stringstream ssrepr;
  std::stringstream ssstr;
  ssrepr << "EnrichedElement(";
  ssstr << "<";

  FiniteElementBaseList::const_iterator it = sub_elements_.begin();
  uint value_size_sum = (*it)->value_shape().size();
  ssrepr << (*it)->repr();
  ssstr << (*it)->str();
  for ( ++it ; it != sub_elements_.end(); ++it)
  {
    ssrepr << ", " << (*it)->repr();
    ssstr << " + " << (*it)->str();

  }

  ssrepr << ")";
  ssstr << ">";
  repr_ = ssrepr.str();
  str_ = ssstr.str();
}

//-----------------------------------------------------------------------------
EnrichedElement::~EnrichedElement()
{
}

//-----------------------------------------------------------------------------
bool const EnrichedElement::is_cellwise_constant() const
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
std::map<uint, uint> const EnrichedElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const EnrichedElement::extract_subelement_component(
    ValueArray const& i) const
{
  return std::pair<uint, uint>();
}

//-----------------------------------------------------------------------------
std::pair<uint, FiniteElementBase const * const> const EnrichedElement::extract_component(ValueArray const& i) const
{
  return std::pair<uint, FiniteElementBase const * const>( i[0] , sub_elements_[i[0]] );
}

//-----------------------------------------------------------------------------
uint const EnrichedElement::num_sub_elements() const
{
  return sub_elements_.size();
}

//-----------------------------------------------------------------------------
FiniteElementBase::FiniteElementBaseList const& EnrichedElement::sub_elements() const
{
  return sub_elements_;
}

}

