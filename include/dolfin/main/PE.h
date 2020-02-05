// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MAIN_PE_H
#define __DOLFIN_MAIN_PE_H

#include <dolfin/main/SubSystemsManager.h>

namespace dolfin
{

/**
 *  @class  PE
 *
 *  @brief  Provides a shim layer for process execution context.
 *
 */

struct PE
{
  /// Returns whether the context is parallel.
  static inline bool parallel()
  {
#if HAVE_MPI
    return SubSystemsManager::active(SubSystemsManager::mpi);
#endif
    return false;
  }

  /// Returns the rank of the calling process in the execution context.
  static uint rank();

  /// Returns the number of processes in the execution context.
  static uint size();

};

} /* namespace dolfin */

#endif /* __DOLFIN_MAIN_PE_H */
