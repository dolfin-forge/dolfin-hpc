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

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

//-----------------------------------------------------------------------------
void dolfin::dolfin_init(int argc, char * argv[])
{
  //--- Process arguments
  uint const maxopt = 1;
  uint curopt = 0;
  int n = 1;
  char const * const * roargv = argv;
  for (int i = 0; i < argc; ++i)
  {
    if(argv[i] == NULL)
    {
      break;
    }
    if(strlen(argv[i]) == 2 && roargv[i][0] == '-')
    {
      int c = argv[i][1];
      char const * argi = roargv[++i];
      switch (c)
        {
        case 'n':
          n = atoi(argi);
          ++curopt;
          break;
        default:
          break;
        }
    }
    if(curopt == maxopt)
    {
      break;
    }
  }

  //--- Initialize subsystems

#ifdef HAVE_MPI
  SubSystemsManager::initMPI(argc, argv, n);
#endif

  // Cannot use MPI functions before initComm
  if (MPI::processNumber() > 0)
  {
    dolfin::LogManager::logger().silence();
  }
  if (MPI::processGlobalNumber() == 0)
  {
    message("Initializing DOLFIN version %s : running on %d %s.\n",
            DOLFIN_VERSION, dolfin::MPI::numGlobalProcesses(),
            (dolfin::MPI::numGlobalProcesses() > 1 ? "processes" : "process"));
  }
  if (MPI::numGroups() > 1)
  {
    message("Group %d/%d : %d %s", dolfin::MPI::groupNumber() + 1,
            dolfin::MPI::numGroups(), dolfin::MPI::numProcesses(),
            (dolfin::MPI::numProcesses() > 1 ? "processes" : "process"));
    if (MPI::groupNumber() > 0)
    {
      dolfin::LogManager::logger().silence();
    }
  }

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
