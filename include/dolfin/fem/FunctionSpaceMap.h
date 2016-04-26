// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2015-03-25
// Last changed: 2015-03-25

#ifndef __DOLFIN_FUNCTION_SPACE_MAP_H
#define __DOLFIN_FUNCTION_SPACE_MAP_H

#include <dolfin/common/types.h>
#include <dolfin/ufl/UFLrepr.h>

#include <map>
#include <string>

namespace ufl
{

  class FiniteElementSpace;

}

namespace dolfin
{

/**
 *  @class  FunctionSpaceMap
 *
 *  @brief  This class serves as provider of function space signature to allow
 *          specifying the type of discrete space to instantiate from generated
 *          source code.
 */

class FunctionSpaceMap
{
  typedef std::map<std::string, ufl::repr> Container;
  typedef std::pair<std::string, ufl::repr> Item;

public:

  ///
  FunctionSpaceMap();

  ///
  ~FunctionSpaceMap();

  ///
  bool has(std::string const& label) const;

  ///
  ufl::repr get(std::string const& label) const;

  ///
  uint size() const;

  ///
  void add(std::string const& label, ufl::FiniteElementSpace const& space);

  ///
  void clear();

private:

  Container spaces_;

};

}

#endif /* __DOLFIN_FUNCTION_SPACE_MAP_H */
