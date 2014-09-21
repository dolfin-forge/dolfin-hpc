// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurélien Larcher, 2014
//
// First added:  2007-04-12
// Last changed: 2014-09-21

#ifndef __ELEMENT_LIBRARY_H
#define __ELEMENT_LIBRARY_H

#include <dolfin/ufl/UFLElementList.h>

#include <ufc.h>

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
