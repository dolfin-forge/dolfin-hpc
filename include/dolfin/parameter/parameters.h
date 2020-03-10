// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PARAMETERS_H
#define __DOLFIN_PARAMETERS_H

#include <dolfin/parameter/Parameter.h>

#include <dolfin/parameter/ParameterSystem.h>

namespace dolfin
{

/// Set value of parameter
void dolfin_set( std::string key, Parameter value );

/// Add parameter
void dolfin_add( std::string key, Parameter value );

/// Get value of parameter with given key
template < typename T >
T dolfin_get( std::string key )
{
	return static_cast< T >( ParameterSystem::parameters.get( key ) );
}

} // end namespace dolfin

#endif
