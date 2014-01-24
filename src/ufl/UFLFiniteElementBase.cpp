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
                                           uint const degree) :
    UFLClass(),
    family_(family),
    cell_(cell),
    degree_(degree),
    value_shape_()
{

  // Check whether cell is valid

  // Check whether degree is valid
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
std::vector<uint> const UFLFiniteElementBase::value_shape() const
{
  return value_shape_;
}

}

