// Copyright (C) 2019 Julian Hornich.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MPI_DATATYPES_H
#define __DOLFIN_MPI_DATATYPES_H

#include <dolfin/common/types.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace dolfin
{

#if !defined( HAVE_MPI )

//-----------------------------------------------------------------------------

typedef void * MPI_Datatype;

template < typename T >
struct MPI_type
{
	static MPI_Datatype const value;
};

template< typename T >
MPI_Datatype const MPI_type<T>::value = nullptr;

#else

//-----------------------------------------------------------------------------

/**
 * @brief Get the MPI Datatype associated to the C/C++ datatype
 *
 * @tparam  T C/C++ datatype
 * @param value MPI_Datatype associated with \c T
 */
template < typename T >
struct MPI_type
{
	// FIXME with c++17 this could be an inline variable, making this a bit cleaner
	static MPI_Datatype const value;
};

//-----------------------------------------------------------------------------

template <>
MPI_Datatype const MPI_type< bool >::value;

// signed integral types
template <>
MPI_Datatype const MPI_type< char >::value;

template <>
MPI_Datatype const MPI_type< signed char >::value;

template <>
MPI_Datatype const MPI_type< signed short int >::value;

template <>
MPI_Datatype const MPI_type< signed int >::value;

template <>
MPI_Datatype const MPI_type< signed long int >::value;

template <>
MPI_Datatype const MPI_type< signed long long >::value;

// unsigned integral types
template <>
MPI_Datatype const MPI_type< unsigned char >::value;

template <>
MPI_Datatype const MPI_type< unsigned short int >::value;

template <>
MPI_Datatype const MPI_type< unsigned long int >::value;

#if ( MPI_VERSION > 1 )
template <>
MPI_Datatype const MPI_type< unsigned long long >::value;
#endif

template <>
MPI_Datatype const MPI_type< unsigned int >::value;

// floating point types
template <>
MPI_Datatype const MPI_type< float >::value;

template <>
MPI_Datatype const MPI_type< double >::value;

template <>
MPI_Datatype const MPI_type< long double >::value;

//-----------------------------------------------------------------------------

#endif

}

#endif
