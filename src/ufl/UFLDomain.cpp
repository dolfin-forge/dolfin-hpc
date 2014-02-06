// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLDomain.h>

#include <dolfin/common/types.h>

#include <iomanip>
#include <iostream>

using dolfin::DOLFIN_UINT_UNDEF;

namespace ufl
{

//-----------------------------------------------------------------------------
Domain::MappingList const Domain::__init_mapping()
{
  MappingList m;
  m.insert(MappingItem(Object::repr_t("None"), Domain::None));
  m.insert(MappingItem(Object::repr_t("'cell1D'"), Domain::cell1D));
  m.insert(MappingItem(Object::repr_t("'cell2D'"), Domain::cell2D));
  m.insert(MappingItem(Object::repr_t("'cell3D'"), Domain::cell3D));
  m.insert(MappingItem(Object::repr_t("'vertex'"), Domain::vertex));
  m.insert(MappingItem(Object::repr_t("'interval'"), Domain::interval));
  m.insert(MappingItem(Object::repr_t("'triangle'"), Domain::triangle));
  m.insert(MappingItem(Object::repr_t("'tetrahedron'"), Domain::tetrahedron));
  m.insert(MappingItem(Object::repr_t("'quadrilateral'"), Domain::quadrilateral));
  m.insert(MappingItem(Object::repr_t("'hexahedron'"), Domain::hexahedron));
  return m;
}

//-----------------------------------------------------------------------------
Domain::DefinitionList const Domain::__init_definitions()
{
  DefinitionList m;
  m.insert(
      DefinitionItem(
          None,
          Definition(DOLFIN_UINT_UNDEF, None, DOLFIN_UINT_UNDEF, "None")));
  m.insert(
      DefinitionItem(cell1D,
                     Definition(1, vertex, DOLFIN_UINT_UNDEF, "'cell1D'")));
  m.insert(
      DefinitionItem(cell2D,
                     Definition(2, cell1D, DOLFIN_UINT_UNDEF, "'cell2D'")));
  m.insert(
      DefinitionItem(cell3D,
                     Definition(3, cell2D, DOLFIN_UINT_UNDEF, "'cell3D'")));
  m.insert(DefinitionItem(vertex, Definition(0, None, 0, "'vertex'")));
  m.insert(DefinitionItem(interval, Definition(1, vertex, 2, "'interval'")));
  m.insert(DefinitionItem(triangle, Definition(2, interval, 3, "'triangle'")));
  m.insert(
      DefinitionItem(tetrahedron, Definition(3, triangle, 4, "'tetrahedron'")));
  m.insert(
      DefinitionItem(quadrilateral,
                     Definition(2, interval, 4, "'quadrilateral'")));
  m.insert(
      DefinitionItem(hexahedron,
                     Definition(3, quadrilateral, 6, "'hexahedron'")));
  return m;
}

//-----------------------------------------------------------------------------
Domain::Domain(Domain::Type const& t) :
    ufl::type<std::string>(Domain::type_repr(t)),
    type_(t)
{
}

//-----------------------------------------------------------------------------
Domain::Domain(repr_t const& repr) :
    ufl::type<std::string>(repr),
    type_(Domain::repr_type(repr))
{
}

//-----------------------------------------------------------------------------
Domain::~Domain()
{
}

//-----------------------------------------------------------------------------
Domain::Type const Domain::type_facet(Domain::Type const& t)
{
  return Definitions().find(t)->second.facet;
}

//-----------------------------------------------------------------------------
uint const Domain::type_dim(Domain::Type const& t)
{
  return Definitions().find(t)->second.dim;
}

//-----------------------------------------------------------------------------
dolfin::uint const Domain::num_facets(Type const& t)
{
  return Definitions().find(t)->second.num_facets;
}

//-----------------------------------------------------------------------------
std::string const Domain::type_repr(Domain::Type const& t)
{
  return Definitions().find(t)->second.str;
}

//-----------------------------------------------------------------------------
Domain::Type const Domain::repr_type(repr_t const& repr)
{
  return Mapping().find(repr)->second;
}

//-----------------------------------------------------------------------------
Domain::Type const Domain::facet() const
{
  return Domain::type_facet(type_);
}

//-----------------------------------------------------------------------------
dolfin::uint const Domain::dim() const
{
  return Domain::type_dim(type_);
}

//-----------------------------------------------------------------------------
dolfin::uint const Domain::num_facets() const
{
  return Domain::type_num_facets(type_);
}

//-----------------------------------------------------------------------------
Domain::Type const Domain::type() const
{
  return type_;
}

//-----------------------------------------------------------------------------
bool const Domain::is_undefined() const
{
  return (type_ == Domain::None);
}

//-----------------------------------------------------------------------------
void Domain::display() const
{
  ufl::type<std::string>::display();
  std::cout << std::setw(24) << "dimension" << " = " << this->dim()
      << std::endl;
  std::cout << std::setw(24) << "facet" << " = "
      << Domain(this->facet()).str() << std::endl;
  std::cout << std::setw(24) << "num_facets" << " = " << this->num_facets()
      << std::endl;
  std::cout << std::endl;
}

}
