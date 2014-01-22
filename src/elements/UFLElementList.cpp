// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/elements/UFLElementList.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
void UFLElementList::register_element(ElementFamilyList& m, FamilyType family,
                                      std::string name, std::string short_name,
                                      uint value_rank, uint degree_min,
                                      uint degree_max, CellTypeSet domains)
{
  ElementDefinition a(name, short_name, value_rank, degree_min, degree_max,
                      domains);
  m.insert(ElementFamilyItem(family, a));
}

//-----------------------------------------------------------------------------
UFLElementList::ElementFamilyList const UFLElementList::__init_elements()
{
  CellTypeSet interval;
  interval.insert(CellType::interval);

  CellTypeSet triangle;
  triangle.insert(CellType::triangle);

  CellTypeSet triangle_tetrahedron;
  triangle_tetrahedron.insert(CellType::triangle);
  triangle_tetrahedron.insert(CellType::tetrahedron);

  CellTypeSet interval_triangle_tetrahedron;
  interval_triangle_tetrahedron.insert(CellType::interval);
  interval_triangle_tetrahedron.insert(CellType::triangle);
  interval_triangle_tetrahedron.insert(CellType::tetrahedron);

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
std::set<CellType::Type> UFLElementList::domains(FamilyType const type) const
{
  return element_definition(type).domains;
}

//-----------------------------------------------------------------------------
void UFLElementList::display() const
{
  for (ElementFamilyList::const_iterator it = Elements.begin();
      it != Elements.end(); ++it)
  {
    std::cout << std::setw(3) << it->first << ":" << std::endl;
    it->second.display();
  }
}

}

