// Copyright (C) 2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_INIT_H
#define __DOLFIN_INIT_H

#include <cstddef>

namespace dolfin
{

/// Initialize DOLFIN (and all configured SubSystems) with command-line
/// arguments. This call is necessary for all
void dolfin_init( int    argc            = 0,
                  char * argv[]          = NULL,
                  long   wallclock_limit = 0,
                  int    parallel_groups = 1 );

void dolfin_finalize();

}

#endif
