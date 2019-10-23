// Copyright (C) 2019 Julian Hornich.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_OPENMP_H
#define __DOLFIN_OPENMP_H

/*
 *  This header enables or disables openmp pragmas
 *  whithout causing the compiler to produce errors
 *
 *  Examples:
 *    - to obtain "#pragma omp atomic" use OPENMP_PRAGMA(atomic)
 *    - to obtain "#pragma omp parallel" use OPENMP_PRAGMA(parallel) or
 *      OPENMP_PARALLEL()
 *    - to obtain "#pragma omp parallel for" use OPENMP_PRAGMA(parallel for) or
 *      OPENMP_PARALLEL(for) or OPENMP_PARALLEL_FOR()
 **/

#include <dolfin/config/dolfin_config.h>

//------------------------------------------------------------------------------
#if defined( HAVE_OPENMP )

#include "omp.h"

// results in #pragma omp <args>
#define OPENMP_PRAGMA( args ) INSERT_PRAGMA( omp args )

// results in #pragma omp parallel <args>
#define OPENMP_PARALLEL( args ) INSERT_PRAGMA( omp parallel args )

// results in #pragma omp parallel for <args>
#define OPENMP_PARALLEL_FOR( args ) INSERT_PRAGMA(omp parallel for args)

#define INSERT_PRAGMA( args ) _Pragma( #args )

#define OPENMP_SCHEDULE schedule( static )

//------------------------------------------------------------------------------
#else

// make void openmp calls and provide dummy functions

#define OPENMP_PRAGMA( args )
#define OPENMP_PARALLEL( args )
#define OPENMP_PARALLEL_FOR( args )
#define OPENMP_SCHEDULE

inline int omp_get_max_threads()
{
	return 1;
}

inline int omp_get_num_threads()
{
	return 1;
}

inline int omp_get_thread_num()
{
	return 0;
}

using omp_lock_t = void *;

inline void omp_init_lock( omp_lock_t t )
{
	( void ) t;
}

inline void omp_set_lock( omp_lock_t t )
{
	( void ) t;
}

inline void omp_unset_lock( omp_lock_t t )
{
	( void ) t;
}

inline void omp_destroy_lock( omp_lock_t t )
{
	( void ) t;
}

inline int omp_test_lock( omp_lock_t t )
{
	( void ) t;
	return 0;
}

inline void omp_set_nested( bool b )
{
	( void ) b;
}

inline bool omp_in_parallel()
{
	return false;
}

#endif

#endif // __DOLFIN_OPENMP_H
