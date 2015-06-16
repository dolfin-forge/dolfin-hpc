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
    List const& elements ) :
    FiniteElementBase("EnrichedElement"),
    sub_elements_(elements),
    family_(Family::Enriched),
    cell_(get_cell(elements)),
    degree_(get_degree_max(elements))
{
  // Create string representation
  std::stringstream ssrepr;
  std::stringstream ssstr;
  ssrepr << "EnrichedElement(";
  ssstr << "<";

  List::const_iterator it = sub_elements_.begin();
  dolfin::uint value_size_sum = (*it)->value_shape().size();
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
Family const& EnrichedElement::family() const
{
  return family_;
}

//-----------------------------------------------------------------------------
Family::Type EnrichedElement::metatype() const
{
  return Family::Enriched;
}

//-----------------------------------------------------------------------------
Cell const& EnrichedElement::cell() const
{
  return cell_;
}

//-----------------------------------------------------------------------------
type<dolfin::uint> const& EnrichedElement::degree() const
{
  return degree_;
}

//-----------------------------------------------------------------------------
ValueArray const& EnrichedElement::value_shape() const
{
  return value_shape_;
}

//-----------------------------------------------------------------------------
bool const EnrichedElement::is_cellwise_constant() const
{
  bool ret = true;
  for ( List::const_iterator it = sub_elements_.begin();
        it != sub_elements_.end(); ++it )
  {
    ret |= (*it)->is_cellwise_constant();
  }
  return ret;
}

//-----------------------------------------------------------------------------
std::map<dolfin::uint, dolfin::uint> const EnrichedElement::symmetry() const
{
  return symmetry_;
}

//-----------------------------------------------------------------------------
std::pair<ValueArray, ValueArray> const EnrichedElement::extract_subelement_component(
    ValueArray const& i) const
{
  return std::pair<ValueArray, ValueArray>();
}

//-----------------------------------------------------------------------------
std::pair<dolfin::uint, FiniteElementBase const *> const EnrichedElement::extract_component(ValueArray const& i) const
{
  return std::pair<dolfin::uint, FiniteElementBase const *>( i[0] , sub_elements_[i[0]] );
}

//-----------------------------------------------------------------------------
dolfin::uint const EnrichedElement::num_sub_elements() const
{
  return sub_elements_.size();
}

//-----------------------------------------------------------------------------
FiniteElementBase::List const& EnrichedElement::sub_elements() const
{
  return sub_elements_;
}

//-----------------------------------------------------------------------------
Object::repr_t const& EnrichedElement::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const& EnrichedElement::str() const
{
  return str_;
}

}

