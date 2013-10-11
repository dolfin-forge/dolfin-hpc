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
  static ufc::finite_element* create_finite_element(std::string const signature)
  {
    return create_finite_element(signature.c_str());
  }

  /// Create dof map with given signature
  static ufc::dof_map* create_dof_map(const char* signature);

  /// Create dof map with given signature
  static ufc::dof_map* create_dof_map(std::string const signature)
  {
    return create_dof_map(signature.c_str());
  }

#if ENABLE_UFL

  static FE::attributes const get_attributes(const char* signature)
  {
    return get_attributes(std::string(signature));
  }

  static FE::attributes const get_attributes(std::string const signature)
  {
    ElementsTable::const_iterator it = Elements.find(signature);
    std::string type;
    std::string family;
    std::string strshape;
    ufc::shape shape;
    uint space = 1;
    uint degree = -1;
    uint value = 1;
    if (it == Elements.end())
    {
      dolfin::warning(
          "Extracting attributes from a finite element not registered in ElementLibrary\n"+signature);
      std::string s(signature);
      size_t pos = s.find("(", 0);
      type = s.substr(0, pos);

      //FIXME: Not mixed element aware
      if (type != FE::MIXED_ELEMENT)
      {
        s.erase(0, pos + 1);
        size_t element = 0;
        size_t t0 = 0;
        size_t t1 = 0;
        std::string tok;
        while ((t1 = s.find(",", t0)) != std::string::npos)
        {
          tok = s.substr(t0, t1 - t0);
          if (element == 0)
          {
            family = s.substr(t0 + 1, t1 - t0 - 2);
          }
          if (element == 1)
          {
            strshape = s.substr(t0 + 7, t1 - t0 - 8);
          }
          if (element == 2)
          {
            std::stringstream ss;
            ss << s.substr(t0 + 7, t1 - t0 - 9);
            ss >> space;
          }
          if (element == 3)
          {
            std::stringstream ss;
            ss << s.substr(t0 + 1, t1 - t0 - 1);
            ss >> degree;
          }
          if (element == 4)
          {
            std::stringstream ss;
            ss << s.substr(t0 + 1, t1 - t0 - 1);
            ss >> value;
          }
          ++element;
          t0 = t1 + 1;
        }
        if (strshape == "interval")
        {
          shape = ufc::interval;
        }
        else if (strshape == "triangle")
        {
          shape = ufc::triangle;
        }
        else if (strshape == "tetrahedron")
        {
          shape = ufc::tetrahedron;
        }
        else
        {
          error("Unknown cell type.");
        }
      }
      else
      {
        family = FE::MIXED_ELEMENT;
      }
      return FE::attributes(type.c_str(), family.c_str(), shape, space, degree,
          value);
    }
    return it->second;
  }

  //-----------------------------------------------------------------------------
  /// List of finite element types
  static Array<std::string> const Types;

//-----------------------------------------------------------------------------
/// List of finite element families
  static Array<std::string> const Families;

//-----------------------------------------------------------------------------
/// Table of finite elements
  typedef std::map<std::string, struct FE::attributes> ElementsTable;
  typedef std::pair<std::string, struct FE::attributes> ElementsItem;
  static ElementsTable const Elements;

#endif

};

}

#endif
