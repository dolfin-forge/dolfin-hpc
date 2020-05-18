// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_DOMAIN_H
#define __DOLFIN_UFL_DOMAIN_H

#include <dolfin/ufl/UFLtype.h>

#include <dolfin/common/types.h>

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
  typedef dolfin::_ordered_set<Domain::Type> Set;

  ///
  Domain(Type const& t);

  ///
  explicit Domain(repr_t const& repr);

  ///
  ~Domain();

  ///
  static Type type_facet(Domain::Type const& t);

  ///
  static dolfin::uint type_dim(Domain::Type const& t);

  ///
  static dolfin::uint type_num_facets(Domain::Type const& t);

  /// UFL:
  Type facet() const;

  /// UFL:
  dolfin::uint dim() const;

  /// UFL:
  dolfin::uint num_facets() const;

  ///
  Type type() const;

  ///
  bool is_undefined() const;

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

    Definition(dolfin::uint a_dim, Type a_facet, dolfin::uint a_num_facets,
               std::string a_str) :
        dim(a_dim),
        facet(a_facet),
        num_facets(a_num_facets),
        str(a_str)
    {
    }
  };

  ///
  static std::string type_repr(Domain::Type const& t);
  typedef dolfin::_ordered_map<Type, Definition> DefinitionList;
  typedef std::pair<Type, Definition> DefinitionItem;
  static DefinitionList const Definitions()
  {
    static DefinitionList const Definitions = __init_definitions();
    return Definitions;
  }
  static DefinitionList const __init_definitions();

  ///
  static Domain::Type repr_type(repr_t const& repr);
  typedef dolfin::_ordered_map<repr_t, Domain::Type> MappingList;
  typedef std::pair<repr_t, Domain::Type> MappingItem;
  static MappingList const Mapping()
  {
    static MappingList const Mapping = __init_mapping();
    return Mapping;
  }
  static MappingList const __init_mapping();

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_DOMAIN_H */
