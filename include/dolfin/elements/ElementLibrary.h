// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_ELEMENT_LIBRARY_H
#define __DOLFIN_ELEMENT_LIBRARY_H

#include <dolfin/ufc/ufc.h>
#include <dolfin/ufl/UFLElementList.h>

#include <string>

namespace dolfin
{

/// Library of pre-generated finite elements and dof maps.

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

  //---------------------------------------------------------------------------
  /// Table of finite elements

private:

  static ufl::ElementList const init_elements();

public:

  static ufl::ElementList const& elements()
  {
    static ufl::ElementList const list = init_elements();
    return list;
  }

};

}

#endif
