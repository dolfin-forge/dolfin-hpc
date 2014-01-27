// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLFiniteElementBase.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLFiniteElementBase::UFLFiniteElementBase(UFLElementList::FamilyType family,
                                           UFLCell const& cell,
                                           uint const degree,
                                           UFLQuadratureScheme quad_scheme,
                                           ValueArray value_shape) :
    UFLClass(),
    family_(family),
    cell_(cell),
    degree_(degree),
    quad_scheme_(quad_scheme),
    value_shape_(value_shape)
{
}

//-----------------------------------------------------------------------------
UFLFiniteElementBase::~UFLFiniteElementBase()
{
}

//-----------------------------------------------------------------------------
UFLElementList::FamilyType const UFLFiniteElementBase::family() const
{
  return family_;
}

//-----------------------------------------------------------------------------
UFLCell const UFLFiniteElementBase::cell() const
{
  return cell_;
}

//-----------------------------------------------------------------------------
uint const UFLFiniteElementBase::degree() const
{
  return degree_;
}

//-----------------------------------------------------------------------------
ValueArray const UFLFiniteElementBase::value_shape() const
{
  return value_shape_;
}

//-----------------------------------------------------------------------------
bool UFLFiniteElementBase::component_is_valid(ValueArray const& i) const
{
  uint r = value_shape_.size();
  bool range_ok = true;
  for(size_t idx = 0; idx < value_shape_.size(); ++idx)
  {
    range_ok = range_ok && (i[idx] < value_shape_[idx]);
  }
  return ( i.size() == r && range_ok);
}

//-----------------------------------------------------------------------------
void UFLFiniteElementBase::check_component(ValueArray const& i) const
{
  if(!component_is_valid(i))
  {
    error("Requested component is invalid");
  }
}

//-----------------------------------------------------------------------------
UFLCell const UFLFiniteElementBase::get_cell(
    FiniteElementBaseList const& elements)
{
  FiniteElementBaseList::const_iterator it = elements.begin();
  UFLCell ret = (*it)->cell();
  for (++it ; it != elements.end(); ++it)
  {
    if( ret.repr() != (*it)->cell().repr())
    {
      error("All subelements of mixed element should have the same cell.");
    }
  }
  return ret;
}

//-----------------------------------------------------------------------------
uint const UFLFiniteElementBase::get_degree_max(
    FiniteElementBaseList const& elements)
{
  uint ret = 0;
  for (FiniteElementBaseList::const_iterator it = elements.begin();
       it != elements.end(); ++it)
  {
    ret = std::max((*it)->degree(), ret);
  }
  return ret;
}

}

