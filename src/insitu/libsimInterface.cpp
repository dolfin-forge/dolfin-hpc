// Copyright (C) 2017 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-05-24
// Last changed: 2017-05-24

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/log.h>
#include <dolfin/insitu/libsimInterface.h>

#include <cstdlib>

using namespace dolfin;

#ifdef HAVE_LIBSIM
//-----------------------------------------------------------------------------
void libsimInterface::initBatch()
{

  if (setupEnv() != VISIT_OKAY)
  {
    error("VisIt/libsim environment initialization error");
  }

  if (VisItInitializeRuntime() != VISIT_OKAY)
  {
    error("VisIt/libsim runtime initialization error");
  }
}
//-----------------------------------------------------------------------------
void libsimInterface::initInteractive()
{
  VisItInitializeSocketAndDumpSimFile("dolfin-hpc", "DOLFIN HPC In-situ viz",
				      "/tmp/", NULL, NULL, NULL);
}
//-----------------------------------------------------------------------------
int libsimInterface::setupEnv()
{
  char *env = NULL;

  if (MPI::rank() == 0)
  {    
    env = VisItGetEnvironment();

    if (env == NULL)
    {
      return VISIT_ERROR;
    }
  }
  
  VisItSetParallel(MPI::size() > 1);
  if (MPI::size() > 1) 
  {
    VisItSetParallelRank(MPI::rank());  
#ifdef HAVE_MPI
    MPI_Comm_dup(MPI::DOLFIN_COMM, &comm);
    VisItSetMPICommunicator((void *) &comm);
#endif
  }

  if (VisItSetupEnvironment2(env) != VISIT_OKAY)
  {
    return VISIT_ERROR;
  }

  if (env != NULL)
  {
    free(env);
  }

  return VISIT_OKAY;

}
//-----------------------------------------------------------------------------
void libsimInterface::shutdown()
{
  VisItDisconnect();
}
//-----------------------------------------------------------------------------
#endif
