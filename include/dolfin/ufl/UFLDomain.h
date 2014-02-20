// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_DOMAIN_H_
#define __UFL_DOMAIN_H_

#include <dolfin/ufl/UFLtype.h>

#include <dolfin/common/types.h>

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
  explicit Domain(repr_t const& repr);

  ///
  ~Domain();

  ///
  static Type const type_facet(Domain::Type const& t);

  ///
  static dolfin::uint const type_dim(Domain::Type const& t);

  ///
  static dolfin::uint const type_num_facets(Domain::Type const& t);

  /// UFL:
  Type const facet() const;

  /// UFL:
  dolfin::uint const dim() const;

  /// UFL:
  dolfin::uint const num_facets() const;

  ///
  Type const type() const;

  ///
  bool const is_undefined() const;

  ///
  void display() const;

private:

  Domain::Type const type_;

  //--- STATIC ----------------------------------------------------------------

  struct Definition
  {
    dolfin::uint dim;
    Type facet;
    dolfin::uint num_facets;
    std::string str;

    Definition(dolfin::uint a_dim, Type a_facet, dolfin::uint a_num_facets, std::string a_str) :
        dim(a_dim),
        facet(a_facet),
        num_facets(a_num_facets),
        str(a_str)
    {
    }
  };

  ///
  static std::string const type_repr(Domain::Type const& t);
  typedef std::map<Type, Definition> DefinitionList;
  typedef std::pair<Type, Definition> DefinitionItem;
  static DefinitionList const Definitions()
  {
    static DefinitionList const Definitions = __init_definitions();
    return Definitions;
  }
  static DefinitionList const __init_definitions();

  ///
  static Domain::Type const repr_type(repr_t const& repr);
  typedef std::map<repr_t const, Domain::Type> MappingList;
  typedef std::pair<repr_t const, Domain::Type> MappingItem;
  static MappingList const Mapping()
  {
    static MappingList const Mapping = __init_mapping();
    return Mapping;
  }
  static MappingList const __init_mapping();

};

} /* namespace ufl */
#endif /* __UFL_DOMAIN_H */
