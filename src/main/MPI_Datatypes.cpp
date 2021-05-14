// Copyright (C) 2019 Julian Hornich.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/main/MPI_Datatypes.h>

namespace dolfin
{

#if defined( DOLFIN_HAVE_MPI )

//-----------------------------------------------------------------------------

template <>
MPI_Datatype const MPI_type< bool >::value = MPI_BYTE;

// signed integral types
template <>
MPI_Datatype const MPI_type< char >::value = MPI_CHAR;

template <>
MPI_Datatype const MPI_type< signed char >::value = MPI_CHAR;

template <>
MPI_Datatype const MPI_type< signed short int >::value = MPI_SHORT;

template <>
MPI_Datatype const MPI_type< signed int >::value = MPI_INT;

template <>
MPI_Datatype const MPI_type< signed long int >::value = MPI_LONG;

template <>
MPI_Datatype const MPI_type< signed long long >::value = MPI_LONG_LONG;

// unsigned integral types
template <>
MPI_Datatype const MPI_type< unsigned char >::value = MPI_UNSIGNED_CHAR;

template <>
MPI_Datatype const MPI_type< unsigned short int >::value = MPI_UNSIGNED_SHORT;

template <>
MPI_Datatype const MPI_type< unsigned long int >::value = MPI_UNSIGNED_LONG;

#if ( MPI_VERSION > 1 )
template <>
MPI_Datatype const MPI_type< unsigned long long >::value = MPI_UNSIGNED_LONG_LONG;
#endif

template <>
MPI_Datatype const MPI_type< unsigned int >::value = MPI_UNSIGNED;

// floating point types
template <>
MPI_Datatype const MPI_type< float >::value = MPI_FLOAT;

template <>
MPI_Datatype const MPI_type< double >::value = MPI_DOUBLE;

template <>
MPI_Datatype const MPI_type< long double >::value = MPI_LONG_DOUBLE;

//-----------------------------------------------------------------------------

#endif

}
