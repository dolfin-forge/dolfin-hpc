// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLCell.h>

#include <dolfin/log/log.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFLCell::UFLCell(UFLDomain::Type const& domain) :
    UFLClass(),
    domain_(domain),
    space_(UFLDomain::dim(domain)),
    invalid_((domain == UFLDomain::None ? true : false)),
    geometric_dimension_(space_.dimension()),
    topological_dimension_(UFLDomain::dim(domain)),
    cell_surface_area_(*this),
    cell_volume_(*this),
    circumradius_(*this),
    facet_area_(*this),
    facet_normal_(*this),
    x_(*this),
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
    topological_dimension_(UFLDomain::dim(domain)),
    cell_surface_area_(*this),
    cell_volume_(*this),
    circumradius_(*this),
    facet_area_(*this),
    facet_normal_(*this),
    x_(*this),
    repr_("Cell(" + UFLDomain::str(domain_) + ", " + space_.repr() + ")"),
    str_("<" + UFLDomain::str(domain_) + " cell in " + space_.str() + ">")
{
}

//-----------------------------------------------------------------------------
UFLSpatialCoordinate const& UFLCell::x() const
{
  return x_;
}

//-----------------------------------------------------------------------------
UFLFacetNormal const& UFLCell::n() const
{
  return facet_normal_;
}

//-----------------------------------------------------------------------------
UFLCellVolume const& UFLCell::volume() const
{
  return cell_volume_;
}

//-----------------------------------------------------------------------------
UFLCircumradius const& UFLCell::circumradius() const
{
  return circumradius_;
}

//-----------------------------------------------------------------------------
UFLFacetArea const& UFLCell::facet_area() const
{
  return facet_area_;
}

//-----------------------------------------------------------------------------
UFLCellSurfaceArea const& UFLCell::surface_area() const
{
  return cell_surface_area_;
}

//-----------------------------------------------------------------------------
bool const& UFLCell::is_undefined() const
{
  return invalid_;
}

//-----------------------------------------------------------------------------
UFLDomain::Type const UFLCell::domain() const
{
  return domain_;
}

//-----------------------------------------------------------------------------
UFLDomain::Type const UFLCell::facet_domain() const
{
  return UFLDomain::facet(domain_);
}

//-----------------------------------------------------------------------------
uint const& UFLCell::num_facets() const
{
  return UFLDomain::num_facets(domain_);
}

//-----------------------------------------------------------------------------
uint const& UFLCell::geometric_dimension() const
{
  return geometric_dimension_;
}

//-----------------------------------------------------------------------------
uint const& UFLCell::topological_dimension() const
{
  return topological_dimension_;
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
