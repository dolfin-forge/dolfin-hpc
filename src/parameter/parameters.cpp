// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2006.
//
// First added:  2005-12-19
// Last changed: 2008-02-11

#include <dolfin/parameter/parameters.h>

#include <dolfin/parameter/ParameterSystem.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Parameter dolfin_get(std::string key)
{
  return ParameterSystem::parameters.get(key);
}
//-----------------------------------------------------------------------------
void dolfin_set(std::string key, Parameter value)
{
  ParameterSystem::parameters.set(key, value);
}
//-----------------------------------------------------------------------------
void dolfin_add(std::string key, Parameter value)
{
  ParameterSystem::parameters.add(key, value);
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
