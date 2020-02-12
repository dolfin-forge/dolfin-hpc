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

typedef void * MPI_Datatype;

template < typename T >
struct MPI_type
{
	static MPI_Datatype value;
};

template< typename T >
MPI_Datatype MPI_type<T>::value = NULL;

#else

/**
 * @brief Get the MPI Datatype associated to the C/C++ datatype
 *
 * @tparam  T C/C++ datatype
 * @param value MPI_Datatype associated with \c T
 */
template < typename T >
struct MPI_type
{
	static MPI_Datatype value;
};

#endif

}

#endif
