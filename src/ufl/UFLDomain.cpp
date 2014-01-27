// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#include <dolfin/ufl/UFLDomain.h>

#include <dolfin/common/types.h>

using dolfin::DOLFIN_UINT_UNDEF;

namespace ufl
{

//-----------------------------------------------------------------------------
Domain::DefinitionList const Domain::__init_domain_definitions()
{
  DefinitionList m;
  m.insert(DefinitionItem(None , Definition(DOLFIN_UINT_UNDEF, None , DOLFIN_UINT_UNDEF, "None") ));
  m.insert(DefinitionItem(cell1D , Definition(1, vertex, DOLFIN_UINT_UNDEF, "cell1D") ));
  m.insert(DefinitionItem(cell2D , Definition(2, cell1D , DOLFIN_UINT_UNDEF, "cell2D") ));
  m.insert(DefinitionItem(cell3D , Definition(3, cell2D , DOLFIN_UINT_UNDEF, "cell3D") ));
  m.insert(DefinitionItem(vertex , Definition(0, None, 0, "vertex") ));
  m.insert(DefinitionItem(interval , Definition(1, vertex, 2, "interval") ));
  m.insert(DefinitionItem(triangle , Definition(2, interval, 3, "triangle") ));
  m.insert(DefinitionItem(tetrahedron , Definition(3, triangle, 4, "tetrahedron") ));
  m.insert(DefinitionItem(quadrilateral , Definition(2, interval, 4, "quadrilateral") ));
  m.insert(DefinitionItem(hexahedron, Definition(3, quadrilateral, 6, "hexahedron") ));
  return m;
}

//-----------------------------------------------------------------------------
Domain::DefinitionList const Domain::DomainDefinitions =
    __init_domain_definitions();

//-----------------------------------------------------------------------------
Domain::Domain(Type const& t) :
    domain_(t),
    repr_(Domain::str(t)),
    str_(Domain::str(t))
{
}

//-----------------------------------------------------------------------------
Domain::~Domain()
{
}

//-----------------------------------------------------------------------------
Domain::Type const Domain::facet(Type const& t)
{
  return DomainDefinitions.find(t)->second.facet;
}

//-----------------------------------------------------------------------------
uint const Domain::dim(Type const& t)
{
  return DomainDefinitions.find(t)->second.dim;
}

//-----------------------------------------------------------------------------
uint const Domain::num_facets(Type const& t)
{
  return DomainDefinitions.find(t)->second.num_facets;
}

//-----------------------------------------------------------------------------
std::string const Domain::str(Type const& t)
{
  return DomainDefinitions.find(t)->second.str;
}

//-----------------------------------------------------------------------------
std::string const Domain::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
std::string const Domain::str() const
{
  return str_;
}

}
