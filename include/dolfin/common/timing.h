// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TIMING_H
#define __DOLFIN_TIMING_H

#include <dolfin/common/types.h>

namespace dolfin
{

/// Start timing
void tic();

/// Return elapsed CPU time
auto toc() -> real;

/// Return and display elapsed CPU time at given verbose level
auto tocd( uint level = 0 ) -> real;

/// Return current CPU time used by process
auto time() -> real;

}

#endif
