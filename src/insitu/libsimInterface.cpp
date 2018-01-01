// Copyright (C) 2017 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-05-24
// Last changed: 2018-01-01

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

  if (setupCallbacks() != VISIT_OKAY)
  {
    error("VisIt/libsim callbacks initialization error");
  }

  InsituData_.batch_ = true;

}
//-----------------------------------------------------------------------------
void libsimInterface::initInteractive()
{

  if (setupEnv() != VISIT_OKAY)
  {
    error("VisIt/libsim environment initialization error");
  }

  VisItInitializeSocketAndDumpSimFile("dolfin-hpc", "DOLFIN HPC In-situ viz",
				      "/tmp/", NULL, NULL, NULL);

  InsituData_.batch_ = false;

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
int libsimInterface::setupCallbacks()
{
  VisItSetGetMetaData(libsimGetMetaData, &InsituData_);
  VisItSetGetDomainList(libsimGetDomain, &InsituData_);
  VisItSetGetVariable(libsimGetFunction, &InsituData_);
  VisItSetGetMesh(libsimGetMesh, &InsituData_);

  return VISIT_OKAY;
}
//-----------------------------------------------------------------------------
void libsimInterface::shutdown()
{
  error("Not implemented yet");

  VisItDisconnect();
}
//-----------------------------------------------------------------------------
void libsimInterface::batchRender(real t, uint tstep)
{
  
  if (!InsituData_.batch_)
  {
    error("VisIt/libsim not initialized in batch mode");
  }

  InsituData_.t_ = t;
  InsituData_.tstep_ = tstep;

  VisItTimeStepChanged();

  VisItUpdatePlots();
  
  // Execute all insitu pipelines
  for(Array<libsimPipeline *>::iterator it = InsituData_.pipelines_.begin();
      it != InsituData_.pipelines_.end(); it++)
  {
    (*it)->exec(InsituData_.t_, InsituData_.tstep_);
  }
  
}
//-----------------------------------------------------------------------------
void libsimInterface::ctrlLoop(real t, uint tstep, int blocking)
{
  if (InsituData_.batch_)
  {
    error("VisIt/libsim not initialized in interactive mode");
  }

  InsituData_.t_ = t;
  InsituData_.tstep_ = tstep;

  VisItTimeStepChanged();

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

    if (setupCallbacks() != VISIT_OKAY)
    {
      error("VisIt/libsim callbacks initialization error");
    }
    

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
void libsimInterface::batchRender(real t, uint tstep)
{
  error("VisIt/libsim is required for in-situ viz");
}
//-----------------------------------------------------------------------------
void libsimInterface::ctrlLoop(real t, uint tstep, int blocking)
{
  error("VisIt/libsim is required for in-situ viz");
}
//-----------------------------------------------------------------------------
#endif
