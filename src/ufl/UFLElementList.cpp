// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLElementList.h>

namespace ufl
{

using dolfin::error;

//-----------------------------------------------------------------------------
void ElementList::register_element(ElementFamilyList& m, Family::Type family,
                                      std::string name, std::string short_name,
                                      uint value_rank, uint degree_min,
                                      uint degree_max, Domain::Set domains)
{
  ElementDefinition a(name, short_name, value_rank, degree_min, degree_max,
                      domains);
  m.insert(ElementFamilyItem(family, a));
}

//-----------------------------------------------------------------------------
ElementList::ElementFamilyList const ElementList::__init_elements()
{
  Domain::Set interval;
  interval.insert(Domain::interval);

  Domain::Set triangle;
  triangle.insert(Domain::triangle);

  Domain::Set triangle_tetrahedron;
  triangle_tetrahedron.insert(Domain::triangle);
  triangle_tetrahedron.insert(Domain::tetrahedron);

  Domain::Set interval_triangle_tetrahedron;
  interval_triangle_tetrahedron.insert(Domain::interval);
  interval_triangle_tetrahedron.insert(Domain::triangle);
  interval_triangle_tetrahedron.insert(Domain::tetrahedron);

  ElementFamilyList m;

  //--- Standard elements -----------------------------------------------------
  register_element(m, Family::ARG, "'Argyris'", "ARG", 0, 1, None,
                   triangle_tetrahedron);

  register_element(m, Family::AW, "'Arnold-Winther'", "AW", 0, 0, None,
                   triangle);

  register_element(m, Family::BDFM, "'Brezzi-Douglas-Fortin-Marini'",
                   "BDFM", 1, 1, None, triangle_tetrahedron);

  register_element(m, Family::BDM, "'Brezzi-Douglas-Marini'", "BDM", 1, 1,
                   None, triangle_tetrahedron);

  register_element(m, Family::CR, "'Crouzeix-Raviart'", "CR", 0, 1, 1,
                   triangle_tetrahedron);

  register_element(m, Family::DG, "'Discontinuous Lagrange'", "DG", 0, 0,
                   None, interval_triangle_tetrahedron);

  register_element(m, Family::HER, "'Hermite'", "HER", 0, 0, None,
                   triangle_tetrahedron);

  register_element(m, Family::CG, "'Lagrange'", "CG", 0, 1, None,
                   interval_triangle_tetrahedron);

  register_element(m, Family::MTW, "'Mardal-Tai-Winther'", "MTW", 0, 0,
                   None, triangle);

  register_element(m, Family::MOR, "'Morley'", "MOR", 0, 0, None,
                   triangle);

  register_element(m, Family::N1curl, "'Nedelec 1st kind H(curl)'",
                   "N1curl", 1, 1, None, triangle_tetrahedron);

  register_element(m, Family::N2curl, "'Nedelec 2nd kind H(curl)'",
                   "N2curl", 1, 1, None, triangle_tetrahedron);

  register_element(m, Family::RT, "'Raviart-Thomas'", "RT", 1, 1, None,
                   triangle_tetrahedron);

  //--- Special elements ------------------------------------------------------
  register_element(m, Family::BQ, "'Boundary Quadrature'", "BQ", 0, 0,
                   None, interval_triangle_tetrahedron);

  register_element(m, Family::B, "'Bubble'", "B", 0, 2, None,
                   interval_triangle_tetrahedron);

  register_element(m, Family::Q, "'Quadrature'", "Q", 0, 0, None,
                   interval_triangle_tetrahedron);

  register_element(m, Family::R, "'Real'", "R", 0, 0, 0,
                   interval_triangle_tetrahedron);

  register_element(m, Family::U, "'Undefined'", "U", 0, 0, None,
                   interval_triangle_tetrahedron);

  return m;
}

//-----------------------------------------------------------------------------
ElementList::ElementFamilyList const ElementList::Elements =
    __init_elements();

//-----------------------------------------------------------------------------
ElementList::ElementList()
{
}

//-----------------------------------------------------------------------------
ElementList::~ElementList()
{
}

//-----------------------------------------------------------------------------
ElementDefinition const ElementList::element_definition(
    Family::Type const type) const
{
  ElementFamilyList::const_iterator it = Elements.find(type);
  if (it == Elements.end())
  {
    error("Required element type has not been found.");
  }
  return it->second;
}

//-----------------------------------------------------------------------------
bool ElementList::has_family(Family::Type const type) const
{
  return (Elements.find(type) != Elements.end());
}

//-----------------------------------------------------------------------------
bool ElementList::has_family_name(std::string const& name) const
{
  ElementFamilyList::const_iterator it = Elements.begin();
  while (it->second.name != name)
  {
    ++it;
  }
  return it != Elements.end();
}

//-----------------------------------------------------------------------------
bool ElementList::has_valid_domain(Family::Type const type,
                                      Domain::Type domain) const
{
  ElementDefinition const d = Elements.find(type)->second;
  return d.domains.count(domain) > 0;
}

//-----------------------------------------------------------------------------
bool ElementList::has_valid_degree(Family::Type const type,
                                      uint const degree) const
{
  ElementDefinition const d = Elements.find(type)->second;
  return (d.degree_range.first <= degree && d.degree_range.second >= degree);
}


//-----------------------------------------------------------------------------
bool ElementList::has_valid_definition(Family::Type const type,
                                          Domain::Type domain,
                                          uint const degree) const
{
  ElementDefinition const d = Elements.find(type)->second;
  return (d.domains.count(domain) > 0 && d.degree_range.first <= degree
          && d.degree_range.second >= degree);
}

//-----------------------------------------------------------------------------
std::string ElementList::name(Family::Type const type) const
{
  return element_definition(type).name;
}

//-----------------------------------------------------------------------------
std::string ElementList::short_name(Family::Type const type) const
{
  return element_definition(type).short_name;
}

//-----------------------------------------------------------------------------
uint ElementList::value_rank(Family::Type const type) const
{
  return element_definition(type).value_rank;
}

//-----------------------------------------------------------------------------
uint ElementList::degree_min(Family::Type const type) const
{
  return element_definition(type).degree_range.first;
}

//-----------------------------------------------------------------------------
uint ElementList::degree_max(Family::Type const type) const
{
  return element_definition(type).degree_range.second;
}

//-----------------------------------------------------------------------------
std::set<Domain::Type> ElementList::domains(Family::Type const type) const
{
  return element_definition(type).domains;
}

//-----------------------------------------------------------------------------
Object::repr_t const ElementList::repr(Family::Type const type) const
{
  return Object::repr_t(name(type));
}

//-----------------------------------------------------------------------------
void ElementList::display() const
{
  for (ElementFamilyList::const_iterator it = Elements.begin();
      it != Elements.end(); ++it)
  {
    it->second.display();
  }
}

}

