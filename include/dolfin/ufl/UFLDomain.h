// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_DOMAIN_H_
#define __UFL_DOMAIN_H_

#include <dolfin/ufl/UFLClass.h>

#include <map>
#include <string>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLDomain
 *
 *  @brief  Provides a C++ equivalent to ufl_domains from ufl.geometry.
 */

class Domain : public Class
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
  Domain(Type const& t);

  ///
  ~Domain();

  ///
  static Type const facet(Type const& t);

  ///
  static uint const dim(Type const& t);

  ///
  static uint const num_facets(Type const& t);

  ///
  static std::string const str(Type const& t);

  /// __repr__
  std::string const repr() const;

  /// __str__
  std::string const str() const;

private:

  Domain::Type const domain_;
  std::string const repr_;
  std::string const str_;

  //--- STATIC ----------------------------------------------------------------

  struct Definition
  {
    uint dim;
    Type facet;
    uint num_facets;
    std::string str;

    Definition(uint a_dim, Type a_facet, uint a_num_facets, std::string a_str) :
        dim(a_dim),
        facet(a_facet),
        num_facets(a_num_facets),
        str(a_str)
    {
    }
  };

  typedef std::map<Type, Definition> DefinitionList;
  typedef std::pair<Type, Definition> DefinitionItem;
  static DefinitionList const DomainDefinitions;
  static DefinitionList const __init_domain_definitions();

};

} /* namespace ufl */
#endif /* __UFL_DOMAIN_H */
