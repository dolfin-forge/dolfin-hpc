// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLCell.h>

#include <dolfin/log/log.h>

#include <iomanip>

namespace ufl
{

//-----------------------------------------------------------------------------
Cell::Cell(Domain const& domain) :
    Class("Cell"),
    domain_(domain),
    space_(domain_.dim()),
    repr_(*this, domain_, space_),
    str_("<" + domain_.str() + " cell in " + space_.str() + ">"),
    invalid_(domain_.is_undefined()),
    geometric_dimension_(space_.dimension()),
    topological_dimension_(domain_.dim()),
    cell_surface_area_(*this),
    cell_volume_(*this),
    circumradius_(*this),
    facet_area_(*this),
    facet_normal_(*this),
    x_(*this)
{
}

//-----------------------------------------------------------------------------
Cell::Cell(Domain const& domain, Space const& space) :
    Class("Cell"),
    domain_(domain),
    space_(space),
    repr_(*this, domain_, space_),
    str_("<" + domain_.str() + " cell in " + space_.str() + ">"),
    invalid_(domain_.is_undefined()),
    geometric_dimension_(space_.dimension()),
    topological_dimension_(domain_.dim()),
    cell_surface_area_(*this),
    cell_volume_(*this),
    circumradius_(*this),
    facet_area_(*this),
    facet_normal_(*this),
    x_(*this)
{
}

//-----------------------------------------------------------------------------
Cell::Cell(repr_t const& repr) :
    Class("Cell", repr),
    domain_(arg(0)),
    space_(arg(1)),
    repr_(*this, domain_, space_),
    str_("<" + domain_.str() + " cell in " + space_.str() + ">"),
    invalid_(domain_.is_undefined()),
    geometric_dimension_(space_.dimension()),
    topological_dimension_(domain_.dim()),
    cell_surface_area_(*this),
    cell_volume_(*this),
    circumradius_(*this),
    facet_area_(*this),
    facet_normal_(*this),
    x_(*this)
{
}

//-----------------------------------------------------------------------------
Cell::~Cell()
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
Domain const Cell::domain() const
{
  return domain_;
}

//-----------------------------------------------------------------------------
Domain const Cell::facet_domain() const
{
  return Domain(domain_.facet());
}

//-----------------------------------------------------------------------------
dolfin::uint const Cell::num_facets() const
{
  return domain_.num_facets();
}

//-----------------------------------------------------------------------------
dolfin::uint const Cell::geometric_dimension() const
{
  return geometric_dimension_;
}

//-----------------------------------------------------------------------------
dolfin::uint const Cell::topological_dimension() const
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

//-----------------------------------------------------------------------------
void Cell::display() const
{
  Class::display();
  std::cout << std::setw(24) << "is_undefined" << " = " << this->is_undefined() << std::endl;
  std::cout << std::setw(24) << "domain" << " = " << this->domain().str() << std::endl;
  std::cout << std::setw(24) << "facet_domain" << " = " << this->facet_domain().str() << std::endl;
  std::cout << std::setw(24) << "geometric_dimension" << " = " << this->geometric_dimension() << std::endl;
  std::cout << std::setw(24) << "topological_dimension" << " = " << this->topological_dimension() << std::endl;
  std::cout << std::endl;
}

}
