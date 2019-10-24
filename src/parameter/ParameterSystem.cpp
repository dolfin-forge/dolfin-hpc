// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <limits>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/parameter/ParameterSystem.h>

// Initialize the global parameter database
dolfin::ParameterSystem dolfin::ParameterSystem::parameters;

using namespace dolfin;

//-----------------------------------------------------------------------------
ParameterSystem::ParameterSystem() : ParameterList()
{
#include <dolfin/parameter/DefaultParameters.h>
}
//-----------------------------------------------------------------------------
