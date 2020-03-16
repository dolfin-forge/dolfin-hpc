// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PARAMETERS_H
#define __DOLFIN_PARAMETERS_H

#include <dolfin/parameter/ParameterSystem.h>

namespace dolfin
{

/// Set value of parameter
template< typename T >
void dolfin_set( std::string key, T const & value );

/// Add parameter
template< typename T >
void dolfin_add( std::string key, T const & value );

/// Get value of parameter with given key
template < typename T >
T dolfin_get( std::string key );

// ----------------------------------------------------------------------------
// implementation of tempalte functions
// ----------------------------------------------------------------------------

/// Set value of parameter
template< typename T >
void dolfin_set( std::string key, T const & value )
{
	ParameterSystem::parameters.set(key, value);
}

/// Add parameter
template< typename T >
void dolfin_add( std::string key, T const & value )
{
	ParameterSystem::parameters.set( key, value );
}

/// Get value of parameter with given key
template < typename T >
T dolfin_get( std::string key )
{
	return ParameterSystem::parameters.get<T>( key );
}

} // end namespace dolfin

#endif
