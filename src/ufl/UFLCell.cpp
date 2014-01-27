// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLCell.h>

#include <dolfin/log/log.h>

namespace ufl
{

//-----------------------------------------------------------------------------
Cell::Cell(Domain::Type const& domain) :
    Class(),
    domain_(domain),
    space_(Domain::dim(domain)),
    invalid_((domain == Domain::None ? true : false)),
    geometric_dimension_(space_.dimension()),
    topological_dimension_(Domain::dim(domain)),
    cell_surface_area_(*this),
    cell_volume_(*this),
    circumradius_(*this),
    facet_area_(*this),
    facet_normal_(*this),
    x_(*this),
    repr_("Cell(" + Domain::str(domain_) + ", " + space_.repr() + ")"),
    str_("<" + Domain::str(domain_) + " cell in " + space_.str() + ">")
{
}

//-----------------------------------------------------------------------------
Cell::Cell(Domain::Type const& domain, Space const& space) :
    Class(),
    domain_(domain),
    space_(space),
    invalid_((domain == Domain::None ? true : false)),
    geometric_dimension_(space_.dimension()),
    topological_dimension_(Domain::dim(domain)),
    cell_surface_area_(*this),
    cell_volume_(*this),
    circumradius_(*this),
    facet_area_(*this),
    facet_normal_(*this),
    x_(*this),
    repr_("Cell(" + Domain::str(domain_) + ", " + space_.repr() + ")"),
    str_("<" + Domain::str(domain_) + " cell in " + space_.str() + ">")
{
}

//-----------------------------------------------------------------------------
SpatialCoordinate const& Cell::x() const
{
  return x_;
}

//-----------------------------------------------------------------------------
FacetNormal const& Cell::n() const
{
  return facet_normal_;
}

//-----------------------------------------------------------------------------
CellVolume const& Cell::volume() const
{
  return cell_volume_;
}

//-----------------------------------------------------------------------------
Circumradius const& Cell::circumradius() const
{
  return circumradius_;
}

//-----------------------------------------------------------------------------
FacetArea const& Cell::facet_area() const
{
  return facet_area_;
}

//-----------------------------------------------------------------------------
CellSurfaceArea const& Cell::surface_area() const
{
  return cell_surface_area_;
}

//-----------------------------------------------------------------------------
bool const& Cell::is_undefined() const
{
  return invalid_;
}

//-----------------------------------------------------------------------------
Domain::Type const& Cell::domain() const
{
  return domain_;
}

//-----------------------------------------------------------------------------
Domain::Type const& Cell::facet_domain() const
{
  return Domain::facet(domain_);
}

//-----------------------------------------------------------------------------
uint const Cell::num_facets() const
{
  return Domain::num_facets(domain_);
}

//-----------------------------------------------------------------------------
uint const& Cell::geometric_dimension() const
{
  return geometric_dimension_;
}

//-----------------------------------------------------------------------------
uint const& Cell::topological_dimension() const
{
  return topological_dimension_;
}

//-----------------------------------------------------------------------------
Object::repr_t const Cell::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const Cell::str() const
{
  return str_;
}

}
