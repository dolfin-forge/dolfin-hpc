// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef  __UFL_ELEMENT_LIST_H_
#define  __UFL_ELEMENT_LIST_H_

#include <dolfin/common/Array.h>
#include <dolfin/mesh/CellType.h>

#include <climits>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLElementList
 *
 *  @brief  Provides the list of elements supported by UFL and specified in file
 *          ufl.elementlist of UFL version 2.1.1.
 *          Only cell types for simplicial meshes have been retained as
 *          dolfin-hpc does not support quadrilateral/hexahedral meshes.
 *
 */

/// Defines None as the maximum integer possible.
#define None INT_MAX
typedef std::set<CellType::Type> CellTypeSet;

struct ElementDefinition
{
  std::string name;
  std::string short_name;
  uint value_rank;
  std::pair<uint, uint> degree_range;
  CellTypeSet domains;

  ElementDefinition(std::string a_name, std::string a_short_name,
                    uint a_value_rank, uint a_degree_min, uint a_degree_max,
                    CellTypeSet set_of_domains) :
      name(a_name),
      short_name(a_short_name),
      value_rank(a_value_rank),
      degree_range(std::pair<uint, uint>(a_degree_min, a_degree_max)),
      domains(set_of_domains)
  {
  }

  void display() const
  {
    std::stringstream ss;
    uint p = 16;
    ss << std::setw(p) << "name" << ": " << name << std::endl;
    ss << std::setw(p) << "short_name" << ": " << short_name << std::endl;
    ss << std::setw(p) << "value_rank" << ": " << value_rank << std::endl;
    ss << std::setw(p) << "degree_min" << ": " << degree_range.first
        << std::endl;
    ss << std::setw(p) << "degree_max" << ": " << degree_range.second
        << std::endl;
    ss << std::setw(p) << "domains" << ": ";
    for (CellTypeSet::const_iterator it = domains.begin(); it != domains.end();
        ++it)
    {
      ss << CellType::type2string(*it) << " ";
    }
    ss << std::endl;
    message(ss.str());
  }
};

class UFLElementList
{

public:

  enum FamilyType
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
    U
  };

  /// Meyers singleton
  static UFLElementList const& Supported()
  {
    static UFLElementList instance_;
    return instance_;
  }

  ///
  bool has_family(FamilyType const type) const;

  ///
  std::string name(FamilyType const type) const;

  ///
  std::string short_name(FamilyType const type) const;

  ///
  uint value_rank(FamilyType const type) const;

  ///
  uint degree_min(FamilyType const type) const;

  ///
  uint degree_max(FamilyType const type) const;

  ///
  std::set<CellType::Type> domains(FamilyType const type) const;

  ///
  void display() const;

private:

  ///
  UFLElementList();

  ///
  ~UFLElementList();

  /// Family list with attributes
  typedef std::pair<FamilyType, ElementDefinition> ElementFamilyItem;
  typedef std::map<FamilyType, ElementDefinition> ElementFamilyList;

  ///
  ElementDefinition const element_definition(FamilyType const type) const;

  static ElementFamilyList const Elements;

  static void register_element(ElementFamilyList& m, FamilyType family,
                               std::string name, std::string short_name,
                               uint value_rank, uint degree_min,
                               uint degree_max, CellTypeSet domains);

  static ElementFamilyList const __init_elements();

};

}

#endif /* __UFL_ELEMENT_LIST_H_ */
