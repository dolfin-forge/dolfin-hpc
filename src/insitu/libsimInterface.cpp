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
void libsimInterface::batchRender(std::string filename)
{
  VisItTimeStepChanged();

  VisItUpdatePlots();
  
  if (VisItSaveWindow(filename.c_str(), 16384, 16384, 
		      VISIT_IMAGEFORMAT_PNG) != VISIT_OKAY)
  {
    error("VisIt failed to render pipeline");
  }
					      
}
//-----------------------------------------------------------------------------
void libsimInterface::ctrlLoop()
{
  int blocking = 0;

  int visit_state = VisItDetectInput(blocking, -1);

  if (visit_state < 0) 
  {
    error("Badness...");
  }
  else if (visit_state == 0)
  {
    return;
  }
  else if (visit_state == 1)
  {
    runflag = 0;
    if (VisItAttemptToCompleteConnection() != VISIT_OKAY)
    {
      error("VisIt failed to connect!");
    }
    message("VisIt connected!");


    while(1) 
    {
      blocking = 1;
      visit_state = VisItDetectInput(blocking, -1);

      if (!VisItProcessEngineCommand())
      {
	VisItDisconnect();
	return;
      }	
    }

  }
}
//-----------------------------------------------------------------------------
#else
void libsimInterface::initBatch()
{
  error("VisIt/libsim is required for in-situ viz");
}
//-----------------------------------------------------------------------------
void libsimInterface::initInteractive()
{
  error("VisIt/libsim is required for in-situ viz");
}
//-----------------------------------------------------------------------------
void libsimInterface::shutdown() 
{
  error("VisIt/libsim is required for in-situ viz");
}
//-----------------------------------------------------------------------------
void libsimInterface::batchRender(std::string filename)
{
  error("VisIt/libsim is required for in-situ viz");
}
//-----------------------------------------------------------------------------
void libsimInterface::ctrlLoop()
{
  error("VisIt/libsim is required for in-situ viz");
}
//-----------------------------------------------------------------------------
#endif
