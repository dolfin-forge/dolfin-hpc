// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson 2009.
//
// First added:  2005-02-13
// Last changed: 2009-04-22

#include <dolfin/main/init.h>

#include <dolfin/common/constants.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/log.h>
#include <dolfin/log/LogManager.h>
#include <dolfin/main/MPI.h>
#include <dolfin/main/SubSystemsManager.h>

//-----------------------------------------------------------------------------
void dolfin::dolfin_init(int argc, char* argv[])
{
#ifdef HAVE_MPI
  SubSystemsManager::initMPI(argc, argv);
#endif

  // Cannot use MPI functions before initComm
  if(MPI::processNumber() > 0)
  {
    dolfin::LogManager::logger().silence();
  }
  message("Initializing DOLFIN version %s : "
          "running on %d/%d %s in group %d/%d.\n", DOLFIN_VERSION,
          dolfin::MPI::numProcesses(),
          dolfin::MPI::numGlobalProcesses(),
          (dolfin::MPI::numProcesses() > 1 ? "procs" : "proc"),
          dolfin::MPI::groupNumber() + 1, dolfin::MPI::numGroups());

#ifdef HAVE_PETSC
  SubSystemsManager::initPETSc(argc, argv);
#endif

#ifdef HAVE_ZOLTAN
  SubSystemsManager::initZoltan(argc, argv);
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
