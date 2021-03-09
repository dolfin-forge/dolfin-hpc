// Copyright (C) 2017 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/log.h>
#include <dolfin/common/types.h>
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
void libsimInterface::init(Mode mode, bool debug)
{

  if (setupEnv() != VISIT_OKAY)
  {
    error("VisIt/libsim environment initialization error");
  }

  if (debug)
  {
    VisItSetOptions("-debug 5 -clobber_vlogs");
  }

  switch(mode)
  {
  case batch:
    if (initBatch() != VISIT_OKAY)
    {
      error("VisIt/libsim batch mode initialization error");
    }
    break;
  case interactive:
    if (initInteractive() != VISIT_OKAY)
    {
      error("VisIt/libsim interactive mode initialization error");
    }
    break;
  default:
    error("Invalid VisIt/libsim mode");
  }
}
//-----------------------------------------------------------------------------
int libsimInterface::initBatch()
{

  if (VisItInitializeRuntime() != VISIT_OKAY)
  {
    error("VisIt/libsim runtime initialization error");
  }

  if (setupCallbacks() != VISIT_OKAY)
  {
    error("VisIt/libsim callbacks initialization error");
  }

  InsituData_.batch_ = true;

  return VISIT_OKAY;

}
//-----------------------------------------------------------------------------
int libsimInterface::initInteractive()
{

  if (VisItInitializeSocketAndDumpSimFile("sim",
					  "DOLFIN HPC In-situ viz",
					  "/tmp",
					  nullptr, nullptr, nullptr) != VISIT_OKAY)
  {
    error("VisIt/libsim socket initialization error");
  }

  InsituData_.batch_ = false;

  return VISIT_OKAY;

}
//-----------------------------------------------------------------------------
int libsimInterface::setupEnv()
{
  char *env = nullptr;
  const std::string visit_path = dolfin_get<std::string>("VisIt directory");
  VisItSetDirectory((char *) visit_path.c_str());

#ifdef DOLFIN_HAVE_MPI
  VisItSetBroadcastIntFunction(libsimBroadcastInt);
  VisItSetBroadcastStringFunction(libsimBroadcastStr);
#endif
  VisItSetParallel(PE::size() > 1);
  VisItSetParallelRank(PE::rank());

#ifdef DOLFIN_HAVE_MPI
  VisItSetMPICommunicator(&MPI::DOLFIN_COMM);
#endif

  if (PE::rank() == 0)
  {
    env = VisItGetEnvironment();

    if (env == nullptr)
    {
      return VISIT_ERROR;
    }
  }

  if (VisItSetupEnvironment2(env) != VISIT_OKAY)
  {
    return VISIT_ERROR;
  }

  if (env != nullptr)
  {
    free(env);
  }

  return VISIT_OKAY;

}
//-----------------------------------------------------------------------------
int libsimInterface::setupCallbacks()
{
#ifdef DOLFIN_HAVE_MPI
  VisItSetSlaveProcessCallback2(libsimSlaveProcess, &InsituData_);
#endif
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
void libsimInterface::batchRender(real t, size_t tstep)
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
  for(std::vector<libsimPipeline *>::iterator it = InsituData_.pipelines_.begin();
      it != InsituData_.pipelines_.end(); it++)
  {
    (*it)->exec(InsituData_.t_, InsituData_.tstep_);
  }

}
//-----------------------------------------------------------------------------
void libsimInterface::ctrlLoop(real t, size_t tstep, int blocking)
{
  if (InsituData_.batch_)
  {
    error("VisIt/libsim not initialized in interactive mode");
  }

  InsituData_.t_ = t;
  InsituData_.tstep_ = tstep;

  VisItTimeStepChanged();

  VisItUpdatePlots();

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
void libsimInterface::init(Mode, bool)
{
  error("VisIt/libsim is required for in-situ viz");
}
//-----------------------------------------------------------------------------
void libsimInterface::shutdown()
{
  error("VisIt/libsim is required for in-situ viz");
}
//-----------------------------------------------------------------------------
void libsimInterface::batchRender(real, size_t)
{
  error("VisIt/libsim is required for in-situ viz");
}
//-----------------------------------------------------------------------------
void libsimInterface::ctrlLoop(real, size_t, int)
{
  error("VisIt/libsim is required for in-situ viz");
}
//-----------------------------------------------------------------------------
#endif
