// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/parameter/parameters.h>

namespace dolfin
{

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
