// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson 2009.
//
// First added:  2005-02-13
// Last changed: 2009-04-22

#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>
#include "SubSystemsManager.h"
#include "MPI.h"
#include "init.h"

//-----------------------------------------------------------------------------
void dolfin::dolfin_init(int argc, char* argv[])
{
  message("Initializing DOLFIN version %s.", DOLFIN_VERSION);
  
#ifdef HAS_PETSC
  SubSystemsManager::initPETSc(argc, argv);
#endif

#ifdef HAS_MPI
  MPI::initComm();
#endif
}
//-----------------------------------------------------------------------------
void dolfin::dolfin_finalize()
{
  // Finalize subsystems in the correct order
  SubSystemsManager::finalizePETSc();
  SubSystemsManager::finalizeMPI();
}
//-----------------------------------------------------------------------------
