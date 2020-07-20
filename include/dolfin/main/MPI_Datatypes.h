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

using MPI_Datatype = void *;

template < typename T >
struct MPI_type
{
	static constexpr MPI_Datatype value = nullptr;
};

//-----------------------------------------------------------------------------

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
	static MPI_Datatype const value;
};

//-----------------------------------------------------------------------------

extern template MPI_Datatype const MPI_type< bool >::value;

//-----------------------------------------------------------------------------

// signed integral types

extern template MPI_Datatype const MPI_type< char >::value;
extern template MPI_Datatype const MPI_type< signed char >::value;
extern template MPI_Datatype const MPI_type< signed short int >::value;
extern template MPI_Datatype const MPI_type< signed int >::value;
extern template MPI_Datatype const MPI_type< signed long int >::value;
extern template MPI_Datatype const MPI_type< signed long long >::value;

//-----------------------------------------------------------------------------

// unsigned integral types

extern template MPI_Datatype const MPI_type< unsigned char >::value;
extern template MPI_Datatype const MPI_type< unsigned short int >::value;
extern template MPI_Datatype const MPI_type< unsigned int >::value;
extern template MPI_Datatype const MPI_type< unsigned long int >::value;

#if ( MPI_VERSION > 1 )
extern template MPI_Datatype const MPI_type< unsigned long long >::value;
#endif

//-----------------------------------------------------------------------------

// floating point types

extern template MPI_Datatype const MPI_type< float >::value;
extern template MPI_Datatype const MPI_type< double >::value;
extern template MPI_Datatype const MPI_type< long double >::value;

//-----------------------------------------------------------------------------

#endif

}

#endif
