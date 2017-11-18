// Copyright (C) 2017 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-05-24
// Last changed: 2017-10-09

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/log.h>
#include <dolfin/common/Array.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/insitu/libsimPipeline.h>
#include <dolfin/insitu/libsimInterface.h>

#include <algorithm>
#include <cstdlib>

using namespace dolfin;

//--- STATIC ------------------------------------------------------------------

int libsimInterface::runflag = 0;
libsimInterface::libsimData libsimInterface::InsituData_;

#ifdef HAVE_MPI
dolfin::MPI::Communicator libsimInterface::comm;
#endif

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

  VisItSetGetMetaData(libsimGetMetaData, &InsituData_);
  VisItSetGetDomainList(libsimGetDomain, &InsituData_);
  VisItSetGetVariable(libsimGetFunction, &InsituData_);
  VisItSetGetMesh(libsimGetMesh, &InsituData_);

}
//-----------------------------------------------------------------------------
void libsimInterface::initInteractive()
{

  error("Not implemented yet");

  if (setupEnv() != VISIT_OKAY)
  {
    error("VisIt/libsim environment initialization error");
  }

  VisItInitializeSocketAndDumpSimFile("dolfin-hpc", "DOLFIN HPC In-situ viz",
				      "/tmp/", NULL, NULL, NULL);

}
//-----------------------------------------------------------------------------
int libsimInterface::setupEnv()
{
  char *env = NULL;
  const std::string visit_path = dolfin_get("VisIt directory");
  VisItSetDirectory((char *) visit_path.c_str());

  VisItSetParallel(PE::size() > 1);
  VisItSetParallelRank(PE::rank());  

  if (PE::rank() == 0)
  {    
    env = VisItGetEnvironment();

    if (env == NULL)
    {
      return VISIT_ERROR;
    }
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
  error("Not implemented yet");

  VisItDisconnect();
}
//-----------------------------------------------------------------------------
void libsimInterface::batchRender()
{
  VisItTimeStepChanged();

  VisItUpdatePlots();
  
  // Execute all insitu pipelines
  for(Array<libsimPipeline *>::iterator it = InsituData_.pipelines_.begin();
      it != InsituData_.pipelines_.end(); it++)
  {
    (*it)->exec();
  }
  
}
//-----------------------------------------------------------------------------
void libsimInterface::ctrlLoop()
{
  error("Not implemented yet");

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
void libsimInterface::batchRender()
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
