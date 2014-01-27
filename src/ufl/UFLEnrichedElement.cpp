// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLEnrichedElement.h>
#include <dolfin/ufl/UFLElementList.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLEnrichedElement::UFLEnrichedElement(
    FiniteElementBaseList const& elements ) :
    UFLFiniteElementBase(UFLElementList::Enriched, get_cell(elements),
                         get_degree_max(elements)),
    sub_elements_(elements)
{
  // Create string representation
  UFLQuadratureScheme qs = "None";

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
UFLEnrichedElement::~UFLEnrichedElement()
{
}

//-----------------------------------------------------------------------------
bool const UFLEnrichedElement::is_cellwise_constant() const
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
std::map<uint, uint> const UFLEnrichedElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const UFLEnrichedElement::extract_subelement_component(
    ValueArray const& i) const
{
  return std::pair<uint, uint>();
}

//-----------------------------------------------------------------------------
std::pair<uint, UFLFiniteElementBase const * const> const UFLEnrichedElement::extract_component(ValueArray const& i) const
{
  return std::pair<uint, UFLFiniteElementBase const * const>( i[0] , sub_elements_[i[0]] );
}

//-----------------------------------------------------------------------------
uint const UFLEnrichedElement::num_sub_elements() const
{
  return sub_elements_.size();
}

//-----------------------------------------------------------------------------
UFLFiniteElementBase::FiniteElementBaseList const& UFLEnrichedElement::sub_elements() const
{
  return sub_elements_;
}

}

