// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-28
// Last changed: 2014-01-28

#ifndef __UFL_FAMILY_H_
#define __UFL_FAMILY_H_

#include <dolfin/ufl/UFLtype.h>
#include <dolfin/ufl/UFLDomain.h>

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

#include <iomanip>
#include <map>
#include <set>
#include <sstream>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLFamily
 *
 *  @brief  Provides a C++ equivalent to family types.
 */

class Family : public type<std::string>
{
  friend class ElementList;

public:

  enum Type
  {
    ARG,
    AW,
    BDFM,
    BDM,
    CR,
    DG,
    HER,
    CG,
    MTW,
    MOR,
    N1curl,
    N2curl,
    RT,
    BQ,
    B,
    Q,
    R,
    U,
    Mixed,
    Vector,
    Tensor,
    Enriched,
    Restricted
  };

  ///
  Family(Family::Type const& t);

  ///
  Family(repr_t const& repr);

  ///
  ~Family();

  ///
  Family::Type const type() const;

  ///
  std::string name() const;

  ///
  std::string short_name() const;

  ///
  dolfin::uint value_rank() const;

  ///
  dolfin::uint degree_min() const;

  ///
  dolfin::uint degree_max() const;

  ///
  Domain::Set domains() const;

  ///
  static bool has_type(Family::Type const type);

  ///
  static bool has_name(std::string const& name);

  ///
  bool has_valid_domain(Domain::Type domain) const;

  ///
  bool has_valid_degree(dolfin::uint const degree) const;

  ///
  bool has_valid_definition(Domain::Type domain, dolfin::uint const degree) const;

private:

  Type const type_;

  struct Definition
  {
    std::string name;
    std::string short_name;
    dolfin::uint value_rank;
    std::pair<dolfin::uint, dolfin::uint> degree_range;
    Domain::Set domains;

    Definition(std::string a_name, std::string a_short_name, dolfin::uint a_value_rank,
               dolfin::uint a_degree_min, dolfin::uint a_degree_max, Domain::Set set_of_domains) :
        name(a_name),
        short_name(a_short_name),
        value_rank(a_value_rank),
        degree_range(std::pair<dolfin::uint, dolfin::uint>(a_degree_min, a_degree_max)),
        domains(set_of_domains)
    {
    }

    void display() const
    {
      std::stringstream ss;
      dolfin::uint p = 16;
      ss << std::setw(p) << "name" << ": " << name << std::endl;
      ss << std::setw(p) << "short_name" << ": " << short_name << std::endl;
      ss << std::setw(p) << "value_rank" << ": " << value_rank << std::endl;
      ss << std::setw(p) << "degree_min" << ": " << degree_range.first
          << std::endl;
      ss << std::setw(p) << "degree_max" << ": " << degree_range.second
          << std::endl;
      ss << std::setw(p) << "domains" << ": ";
      for (Domain::Set::const_iterator it = domains.begin();
          it != domains.end(); ++it)
      {
        ss << Domain(*it).str() << " ";
      }
      ss << std::endl;
      dolfin::message(ss.str());
    }
  };

  Definition const def_;

  ///
  static std::string const type_repr(Family::Type const& t);
  typedef std::map<Type, Definition> DefinitionList;
  typedef std::pair<Type, Definition> DefinitionItem;
  static dolfin::uint const None = dolfin::DOLFIN_UINT_UNDEF;
  static DefinitionList const Definitions()
  {
    static DefinitionList const definitions = __init_definitions();
    return definitions;
  }
  static DefinitionList const __init_definitions();

  static void register_family(DefinitionList& m, Family::Type family,
                              std::string name, std::string short_name,
                              dolfin::uint value_rank, dolfin::uint degree_min, dolfin::uint degree_max,
                              Domain::Set domains);

  ///
  static Family::Type const repr_type(repr_t const& repr);
#ifdef __SUNPRO_CC
  typedef std::map<repr_t, Family::Type> MappingList;
  typedef std::pair<repr_t, Family::Type> MappingItem;
#else
  typedef std::map<repr_t const, Family::Type> MappingList;
  typedef std::pair<repr_t const, Family::Type> MappingItem;
#endif
  static MappingList const Mapping()
  {
    static MappingList const mapping = __init_mapping();
    return mapping;
  }
  static MappingList const __init_mapping();

  ///
  Family::Definition const element_definition(Family::Type const type) const;

};

} /* namespace ufl */
#endif /* __UFL_FAMILY_H_ */
