// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/ufl/UFLFamily.h>

#include <iomanip>
#include <iostream>

namespace ufl
{

//-----------------------------------------------------------------------------
Family::MappingList Family::__init_mapping()
{
  MappingList m;
  m.insert(MappingItem(Object::repr_t("'Argyris'"), Family::ARG));
  m.insert(MappingItem(Object::repr_t("'Arnold-Winther'"), Family::AW));
  m.insert(
      MappingItem(Object::repr_t("'Brezzi-Douglas-Fortin-Marini'"),
                  Family::BDFM));
  m.insert(MappingItem(Object::repr_t("'Brezzi-Douglas-Marini'"), Family::BDM));
  m.insert(MappingItem(Object::repr_t("'Crouzeix-Raviart'"), Family::CR));
  m.insert(MappingItem(Object::repr_t("'Discontinuous Lagrange'"), Family::DG));
  m.insert(MappingItem(Object::repr_t("'Hermite'"), Family::HER));
  m.insert(MappingItem(Object::repr_t("'Lagrange'"), Family::CG));
  m.insert(MappingItem(Object::repr_t("'Mardal-Tai-Winther'"), Family::MTW));
  m.insert(MappingItem(Object::repr_t("'Morley'"), Family::MOR));
  m.insert(
      MappingItem(Object::repr_t("'Nedelec 1st kind H(curl)'"),
                  Family::N1curl));
  m.insert(
      MappingItem(Object::repr_t("'Nedelec 2nd kind H(curl)'"),
                  Family::N2curl));
  m.insert(MappingItem(Object::repr_t("'Raviart-Thomas'"), Family::RT));
  m.insert(MappingItem(Object::repr_t("'Boundary Quadrature'"), Family::BQ));
  m.insert(MappingItem(Object::repr_t("'Bubble'"), Family::B));
  m.insert(MappingItem(Object::repr_t("'Quadrature'"), Family::Q));
  m.insert(MappingItem(Object::repr_t("'Real'"), Family::R));
  m.insert(MappingItem(Object::repr_t("'Undefined'"), Family::U));
  m.insert(MappingItem(Object::repr_t("'Mixed'"), Family::Mixed));
  m.insert(MappingItem(Object::repr_t("'Vector'"), Family::Vector));
  m.insert(MappingItem(Object::repr_t("'Tensor'"), Family::Tensor));
  m.insert(MappingItem(Object::repr_t("'Enriched'"), Family::Enriched));
  m.insert(MappingItem(Object::repr_t("'Restricted'"), Family::Restricted));
  return m;
}

//-----------------------------------------------------------------------------
void Family::register_family(DefinitionList& m, Family::Type family,
                             std::string name, std::string short_name,
                             dolfin::uint value_rank, dolfin::uint degree_min,
                             dolfin::uint degree_max, Domain::Set domains)
{
  Definition a(name, short_name, value_rank, degree_min, degree_max, domains);
  m.insert(DefinitionItem(family, a));
}

//-----------------------------------------------------------------------------
Family::DefinitionList const Family::__init_definitions()
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

  DefinitionList m;

  //--- Standard elements -----------------------------------------------------
  register_family(m, Family::ARG, "'Argyris'", "ARG", 0, 1, None,
                  triangle_tetrahedron);
  register_family(m, Family::AW, "'Arnold-Winther'", "AW", 0, 0, None,
                  triangle);
  register_family(m, Family::BDFM, "'Brezzi-Douglas-Fortin-Marini'", "BDFM", 1,
                  1, None, triangle_tetrahedron);
  register_family(m, Family::BDM, "'Brezzi-Douglas-Marini'", "BDM", 1, 1, None,
                  triangle_tetrahedron);
  register_family(m, Family::CR, "'Crouzeix-Raviart'", "CR", 0, 1, 1,
                  triangle_tetrahedron);
  register_family(m, Family::DG, "'Discontinuous Lagrange'", "DG", 0, 0, None,
                  interval_triangle_tetrahedron);
  register_family(m, Family::HER, "'Hermite'", "HER", 0, 0, None,
                  triangle_tetrahedron);
  register_family(m, Family::CG, "'Lagrange'", "CG", 0, 1, None,
                  interval_triangle_tetrahedron);
  register_family(m, Family::MTW, "'Mardal-Tai-Winther'", "MTW", 0, 0, None,
                  triangle);
  register_family(m, Family::MOR, "'Morley'", "MOR", 0, 0, None, triangle);
  register_family(m, Family::N1curl, "'Nedelec 1st kind H(curl)'", "N1curl", 1,
                  1, None, triangle_tetrahedron);
  register_family(m, Family::N2curl, "'Nedelec 2nd kind H(curl)'", "N2curl", 1,
                  1, None, triangle_tetrahedron);
  register_family(m, Family::RT, "'Raviart-Thomas'", "RT", 1, 1, None,
                  triangle_tetrahedron);

  //--- Special elements ------------------------------------------------------
  register_family(m, Family::BQ, "'Boundary Quadrature'", "BQ", 0, 0, None,
                  interval_triangle_tetrahedron);
  register_family(m, Family::B, "'Bubble'", "B", 0, 2, None,
                  interval_triangle_tetrahedron);
  register_family(m, Family::Q, "'Quadrature'", "Q", 0, 0, None,
                  interval_triangle_tetrahedron);
  register_family(m, Family::R, "'Real'", "R", 0, 0, 0,
                  interval_triangle_tetrahedron);
  register_family(m, Family::U, "'Undefined'", "U", 0, 0, None,
                  interval_triangle_tetrahedron);

  //--- Meta elements ---------------------------------------------------------
  register_family(m, Family::Mixed, "'Mixed'", "MIXED", 0, 0, None,
                  interval_triangle_tetrahedron);
  register_family(m, Family::Vector, "'Vector'", "VECTOR", 0, 0, None,
                  interval_triangle_tetrahedron);
  register_family(m, Family::Tensor, "'Tensor'", "TENSOR", 0, 0, None,
                  interval_triangle_tetrahedron);
  register_family(m, Family::Enriched, "'Enriched'", "ENRICHED", 0, 0, None,
                  interval_triangle_tetrahedron);
  register_family(m, Family::Restricted, "'Restricted'", "RESTRICTED", 0, 0,
                  None, interval_triangle_tetrahedron);

  return m;
}

//-----------------------------------------------------------------------------
Family::Family(Family::Type const& t) :
    ufl::type<std::string>(Family::type_repr(t)),
    type_(t),
    def_(Definitions().find(type_)->second)
{
}

//-----------------------------------------------------------------------------
Family::Family(repr_t const& repr) :
    ufl::type<std::string>(repr),
    type_(Family::repr_type(repr)),
    def_(Definitions().find(type_)->second)
{
}

//-----------------------------------------------------------------------------
Family::~Family()
{
}

//-----------------------------------------------------------------------------
Family::Type Family::type() const
{
  return type_;
}

//-----------------------------------------------------------------------------
std::string Family::name() const
{
  return def_.name;
}

//-----------------------------------------------------------------------------
std::string Family::short_name() const
{
  return def_.short_name;
}

//-----------------------------------------------------------------------------
dolfin::uint Family::value_rank() const
{
  return def_.value_rank;
}

//-----------------------------------------------------------------------------
dolfin::uint Family::degree_min() const
{
  return def_.degree_range.first;
}

//-----------------------------------------------------------------------------
dolfin::uint Family::degree_max() const
{
  return def_.degree_range.second;
}

//-----------------------------------------------------------------------------
dolfin::_ordered_set<Domain::Type> Family::domains() const
{
  return def_.domains;
}

//-----------------------------------------------------------------------------
Family::Definition Family::element_definition(Family::Type const type) const
{
  Family::DefinitionList::const_iterator it = Family::Definitions().find(type);
  if (it == Family::Definitions().end())
  {
    dolfin::error("Required element type has not been found.");
  }
  return it->second;
}

//-----------------------------------------------------------------------------
bool Family::has_type(Family::Type const type)
{
  return (Family::Definitions().find(type) != Family::Definitions().end());
}

//-----------------------------------------------------------------------------
bool Family::has_name(std::string const& name)
{
  Family::DefinitionList::const_iterator it = Family::Definitions().begin();
  while (it->second.name != name)
  {
    ++it;
  }
  return it != Family::Definitions().end();
}

//-----------------------------------------------------------------------------
bool Family::has_valid_domain(Domain::Type domain) const
{
  return def_.domains.count(domain) > 0;
}

//-----------------------------------------------------------------------------
bool Family::has_valid_degree(dolfin::uint const degree) const
{
  return (def_.degree_range.first <= degree
      && def_.degree_range.second >= degree);
}

//-----------------------------------------------------------------------------
bool Family::has_valid_definition(Domain::Type domain,
                                  dolfin::uint const degree) const
{
  return (has_valid_domain(domain) && has_valid_degree(degree));
}

//-----------------------------------------------------------------------------
std::string Family::type_repr(Family::Type const& t)
{
  return Definitions().find(t)->second.name;
}

//-----------------------------------------------------------------------------
Family::Type Family::repr_type(repr_t const& repr)
{
  return Mapping().find(repr)->second;
}

} // end namespace ufl
