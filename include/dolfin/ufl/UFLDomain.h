// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_DOMAIN_H_
#define __UFL_DOMAIN_H_

#include <map>
#include <string>

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLDomain
 *
 *  @brief  Provides a C++ equivalent to ufl_domains from ufl.geometry.
 */

class UFLDomain
{

public:

  enum Type
  {
    None,
    cell1D,
    cell2D,
    cell3D,
    vertex,
    interval,
    triangle,
    tetrahedron,
    quadrilateral,
    hexahedron
  };

  ///
  static Type const facet(Type const& t);

  ///
  static uint const dim(Type const& t);

  ///
  static uint const num_facet(Type const& t);

  ///
  static std::string const str(Type const& t);

private:

  UFLDomain()
  {
  }

  ~UFLDomain()
  {
  }

  struct Definition
  {
    uint dim;
    Type facet;
    uint num_facet;
    std::string str;

    Definition(uint a_dim, Type a_facet, uint a_num_facet, std::string a_str) :
        dim(a_dim),
        facet(a_facet),
        num_facet(a_num_facet),
        str(a_str)
    {
    }
  };

  typedef std::map<Type, Definition> DefinitionList;
  typedef std::pair<Type, Definition> DefinitionItem;
  static DefinitionList const DomainDefinitions;
  static DefinitionList const __init_domain_definitions();

};

} /* namespace dolfin */
#endif /* __UFL_DOMAIN_H */
