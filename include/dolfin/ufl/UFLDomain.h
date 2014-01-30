// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_DOMAIN_H_
#define __UFL_DOMAIN_H_

#include <dolfin/ufl/UFLtype.h>

#include <map>
#include <set>
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

class Domain : public type<std::string>
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
  typedef std::set<Domain::Type> Set;

  ///
  Domain(Type const& t);

  ///
  ~Domain();

  ///
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

  static DefinitionList const Definitions()
  {
    static DefinitionList const DomainDefinitions = __init_definitions();
    return DomainDefinitions;
  }

  ///
  static Type const type_facet(Domain::Type const& t);

  ///
  static uint const type_dim(Domain::Type const& t);

  ///
  static uint const type_num_facets(Domain::Type const& t);

  ///
  static std::string const type_str(Domain::Type const& t);

  /// UFL:
  Type const facet() const;

  /// UFL:
  uint const dim() const;

  /// UFL:
  uint const num_facets() const;

  ///
  Type const type() const;

  ///
  bool const is_undefined() const;

  ///
  void display() const;

private:

  Domain::Type const type_;

  //--- STATIC ----------------------------------------------------------------
  static DefinitionList const __init_definitions();

};

} /* namespace ufl */
#endif /* __UFL_DOMAIN_H */
