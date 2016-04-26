// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2015-03-25
// Last changed: 2015-03-25

#include <dolfin/fem/FunctionSpaceMap.h>

#include <dolfin/ufl/UFLFiniteElementSpace.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
FunctionSpaceMap::FunctionSpaceMap()
{
}

//-----------------------------------------------------------------------------
FunctionSpaceMap::~FunctionSpaceMap()
{
}

//-----------------------------------------------------------------------------
bool FunctionSpaceMap::has(std::string const& label) const
{
  return (spaces_.count(label) > 0);
}

//-----------------------------------------------------------------------------
ufl::repr FunctionSpaceMap::get(std::string const& label) const
{
  Container::const_iterator it = spaces_.find(label);
  if (spaces_.find(label) == spaces_.end())
  {
    error("Space label '%s' not defined in function space map.", label.c_str());
  }
  return it->second;
}

//-----------------------------------------------------------------------------
uint FunctionSpaceMap::size() const
{
  return spaces_.size();
}

//-----------------------------------------------------------------------------
void FunctionSpaceMap::add(std::string const& label,
                           ufl::FiniteElementSpace const& space)
{
  if (this->has(label))
  {
    error("Space label '%s' already defined in function space map.",
          label.c_str());
  }
  spaces_.insert(Item(label, space.repr()));
}

//-----------------------------------------------------------------------------
void FunctionSpaceMap::clear()
{
  spaces_.clear();
}

}

