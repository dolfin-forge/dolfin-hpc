// Copyright (C) 2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_INIT_H
#define __DOLFIN_INIT_H

#include <dolfin/common/types.h>

namespace dolfin
{
  
  /// Initialize DOLFIN (and PETSc) with command-line arguments. This
  /// should not be needed in most cases since the initialization is
  /// otherwise handled automatically.
  void dolfin_init(int argc = 0, char* argv[] = nullptr, 
		   long w_limit = 0, int n = 1);
  void dolfin_finalize();

}

#endif
