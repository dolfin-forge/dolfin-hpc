// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLCell.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLCell::UFLCell(UFLDomain::Type const& domain) :
    UFLClass(),
    domain_(domain),
    space_(UFLDomain::dim(domain)),
    invalid_((domain == UFLDomain::None ? true : false)),
    geometric_dimension_(space_.dimension()),
    topological_dimension_(),
    repr_("Cell(" + UFLDomain::str(domain_) + ", " + space_.repr() + ")"),
    str_("<" + UFLDomain::str(domain_) + " cell in " + space_.str() + ">")
{
}

//-----------------------------------------------------------------------------
UFLCell::UFLCell(UFLDomain::Type const& domain, UFLSpace const& space) :
    UFLClass(),
    domain_(domain),
    space_(space),
    invalid_((domain == UFLDomain::None ? true : false)),
    geometric_dimension_(space_.dimension()),
    topological_dimension_(),
    repr_("Cell(" + UFLDomain::str(domain_) + ", " + space_.repr() + ")"),
    str_("<" + UFLDomain::str(domain_) + " cell in " + space_.str() + ">")
{
}

//-----------------------------------------------------------------------------
std::string const UFLCell::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const UFLCell::str() const
{
  return str_;
}

}
