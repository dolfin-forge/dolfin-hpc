// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-04-12
// Last changed: 2007-04-13

#ifndef __ELEMENT_LIBRARY_H
#define __ELEMENT_LIBRARY_H

#include <dolfin/common/Array.h>
#include <dolfin/elements/FE.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/CellType.h>
#include <ufc.h>

#include <string>

namespace dolfin
{

/// Library of pregenerated finite elements and dof maps.

class ElementLibrary
{
public:

  /// Create finite element with given signature
  static ufc::finite_element* create_finite_element(const char* signature);

  /// Create finite element with given signature
  static ufc::finite_element* create_finite_element(std::string const signature);

  /// Create dof map with given signature
  static ufc::dofmap* create_dof_map(const char* signature);

  /// Create dof map with given signature
  static ufc::dofmap* create_dof_map(std::string const signature);

#if ENABLE_UFL

  ///
  static FE::attributes const get_attributes(const char* signature);

  ///
  static FE::attributes const get_attributes(std::string const signature);

  //---------------------------------------------------------------------------
  /// Table of finite elements
  typedef std::map<std::string, struct FE::attributes> ElementsTable;
  typedef std::pair<std::string, struct FE::attributes> ElementsItem;
  static ElementsTable const Elements;

#endif

};

}

#endif
