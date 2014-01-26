// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLElementList.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
void UFLElementList::register_element(ElementFamilyList& m, FamilyType family,
                                      std::string name, std::string short_name,
                                      uint value_rank, uint degree_min,
                                      uint degree_max, DomainSet domains)
{
  ElementDefinition a(name, short_name, value_rank, degree_min, degree_max,
                      domains);
  m.insert(ElementFamilyItem(family, a));
}

//-----------------------------------------------------------------------------
UFLElementList::ElementFamilyList const UFLElementList::__init_elements()
{
  DomainSet interval;
  interval.insert(UFLDomain::interval);

  DomainSet triangle;
  triangle.insert(UFLDomain::triangle);

  DomainSet triangle_tetrahedron;
  triangle_tetrahedron.insert(UFLDomain::triangle);
  triangle_tetrahedron.insert(UFLDomain::tetrahedron);

  DomainSet interval_triangle_tetrahedron;
  interval_triangle_tetrahedron.insert(UFLDomain::interval);
  interval_triangle_tetrahedron.insert(UFLDomain::triangle);
  interval_triangle_tetrahedron.insert(UFLDomain::tetrahedron);

  ElementFamilyList m;

  //--- Standard elements -----------------------------------------------------
  register_element(m, UFLElementList::ARG, "Argyris", "ARG", 0, 1, None,
                   triangle_tetrahedron);

  register_element(m, UFLElementList::AW, "Arnold-Winther", "AW", 0, 0, None,
                   triangle);

  register_element(m, UFLElementList::BDFM, "Brezzi-Douglas-Fortin-Marini",
                   "BDFM", 1, 1, None, triangle_tetrahedron);

  register_element(m, UFLElementList::BDM, "Brezzi-Douglas-Marini", "BDM", 1, 1,
                   None, triangle_tetrahedron);

  register_element(m, UFLElementList::CR, "Crouzeix-Raviart", "CR", 0, 1, 1,
                   triangle_tetrahedron);

  register_element(m, UFLElementList::DG, "Discontinuous Lagrange", "DG", 0, 0,
                   None, interval_triangle_tetrahedron);

  register_element(m, UFLElementList::HER, "Hermite", "HER", 0, 0, None,
                   triangle_tetrahedron);

  register_element(m, UFLElementList::CG, "Lagrange", "CG", 0, 1, None,
                   interval_triangle_tetrahedron);

  register_element(m, UFLElementList::MTW, "Mardal-Tai-Winther", "MTW", 0, 0,
                   None, triangle);

  register_element(m, UFLElementList::MOR, "Morley", "MOR", 0, 0, None,
                   triangle);

  register_element(m, UFLElementList::N1curl, "Nedelec 1st kind H(curl)",
                   "N1curl", 1, 1, None, triangle_tetrahedron);

  register_element(m, UFLElementList::N2curl, "Nedelec 2nd kind H(curl)",
                   "N2curl", 1, 1, None, triangle_tetrahedron);

  register_element(m, UFLElementList::RT, "Raviart-Thomas", "RT", 1, 1, None,
                   triangle_tetrahedron);

  //--- Special elements ------------------------------------------------------
  register_element(m, UFLElementList::BQ, "Boundary Quadrature", "BQ", 0, 0,
                   None, interval_triangle_tetrahedron);

  register_element(m, UFLElementList::B, "Bubble", "B", 0, 2, None,
                   interval_triangle_tetrahedron);

  register_element(m, UFLElementList::Q, "Quadrature", "Q", 0, 0, None,
                   interval_triangle_tetrahedron);

  register_element(m, UFLElementList::R, "Real", "R", 0, 0, 0,
                   interval_triangle_tetrahedron);

  register_element(m, UFLElementList::U, "Undefined", "U", 0, 0, None,
                   interval_triangle_tetrahedron);

  return m;
}

//-----------------------------------------------------------------------------
UFLElementList::ElementFamilyList const UFLElementList::Elements =
    __init_elements();

//-----------------------------------------------------------------------------
UFLElementList::UFLElementList()
{
}

//-----------------------------------------------------------------------------
UFLElementList::~UFLElementList()
{
}

//-----------------------------------------------------------------------------
ElementDefinition const UFLElementList::element_definition(
    FamilyType const type) const
{
  ElementFamilyList::const_iterator it = Elements.find(type);
  if (it == Elements.end())
  {
    error("Required element type has not been found.");
  }
  return it->second;
}

//-----------------------------------------------------------------------------
bool UFLElementList::has_family(FamilyType const type) const
{
  return (Elements.find(type) != Elements.end());
}

//-----------------------------------------------------------------------------
bool UFLElementList::has_family_name(std::string const& name) const
{
  ElementFamilyList::const_iterator it = Elements.begin();
  while (it->second.name != name)
  {
    ++it;
  }
  return it != Elements.end();
}

//-----------------------------------------------------------------------------
bool UFLElementList::has_valid_domain(FamilyType const type,
                                      UFLDomain::Type domain) const
{
  ElementDefinition const d = Elements.find(type)->second;
  return d.domains.count(domain) > 0;
}

//-----------------------------------------------------------------------------
bool UFLElementList::has_valid_degree(FamilyType const type,
                                      uint const degree) const
{
  ElementDefinition const d = Elements.find(type)->second;
  return (d.degree_range.first <= degree && d.degree_range.second >= degree);
}


//-----------------------------------------------------------------------------
bool UFLElementList::has_valid_definition(FamilyType const type,
                                          UFLDomain::Type domain,
                                          uint const degree) const
{
  ElementDefinition const d = Elements.find(type)->second;
  return (d.domains.count(domain) > 0 && d.degree_range.first <= degree
          && d.degree_range.second >= degree);
}

//-----------------------------------------------------------------------------
std::string UFLElementList::name(FamilyType const type) const
{
  return element_definition(type).name;
}

//-----------------------------------------------------------------------------
std::string UFLElementList::short_name(FamilyType const type) const
{
  return element_definition(type).short_name;
}

//-----------------------------------------------------------------------------
uint UFLElementList::value_rank(FamilyType const type) const
{
  return element_definition(type).value_rank;
}

//-----------------------------------------------------------------------------
uint UFLElementList::degree_min(FamilyType const type) const
{
  return element_definition(type).degree_range.first;
}

//-----------------------------------------------------------------------------
uint UFLElementList::degree_max(FamilyType const type) const
{
  return element_definition(type).degree_range.second;
}

//-----------------------------------------------------------------------------
std::set<UFLDomain::Type> UFLElementList::domains(FamilyType const type) const
{
  return element_definition(type).domains;
}

//-----------------------------------------------------------------------------
void UFLElementList::display() const
{
  for (ElementFamilyList::const_iterator it = Elements.begin();
      it != Elements.end(); ++it)
  {
    it->second.display();
  }
}

}

