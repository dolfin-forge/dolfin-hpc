// Copyright (C) 2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Niclas Jansson 2009-2015.
//
// First added:  2008-01-07
// Last changed: 2015-01-30

#include <dolfin/common/constants.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/main/SubSystemsManager.h>

#ifdef HAVE_PETSC
#include <petsc.h>
#endif

#ifdef HAVE_SLEPC
#include <slepc.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_ZOLTAN
#include <zoltan_cpp.h>
#endif

namespace dolfin
{

// Initialise static data
SubSystemsManager SubSystemsManager::sub_systems_manager;
uint SubSystemsManager::mpi_init_sema_ = 0;

//-----------------------------------------------------------------------------
SubSystemsManager::SubSystemsManager() :
    petsc_initialized(false),
    petsc_controls_mpi(false),
    zoltan_initialized(false)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SubSystemsManager::SubSystemsManager(SubSystemsManager const& sub_sys_manager)
{
  error("Should not be using copy constructor of SubSystemsManager.");
}
//-----------------------------------------------------------------------------
SubSystemsManager::~SubSystemsManager()
{
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::initMPI(int argc, char* argv[], uint n)
{
  ++mpi_init_sema_;

#ifdef HAVE_MPI
  if (MPIinitialized())
  {
    return false;
  }

#ifdef HAVE_JANPACK_MPI
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
  MPI_Init(&argc, &argv);
#endif

  MPI::initComm(n);

#else
  // Do nothing
#endif

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::initPETSc()
{
#ifdef HAVE_PETSC
  if (sub_systems_manager.petsc_initialized)
  {
    return false;
  }

  message(1, "Initializing PETSc (ignoring command-line arguments).");

  // Dummy command-line arguments for PETSc. This is needed since
  // PetscInitializeNoArguments() does not seem to work.
  int argc = 0;
  char** argv = NULL;

  // Initialize PETSc
  return initPETSc(argc, argv, false);

#else
  error("DOLFIN has not been configured for PETSc.");
  return true;

#endif
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::initPETSc(int argc, char* argv[], bool cmd_line_args)
{
#ifdef HAVE_PETSC
  if (sub_systems_manager.petsc_initialized)
  {
    return false;
  }

  // Get status of MPI before PETSc initialisation
  const bool mpi_init_status = MPIinitialized();

  // FIXME: What does this do?
  if (cmd_line_args)
  {
    message(1, "Initializing PETSc with given command-line arguments.");
  }

  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

#ifdef HAVE_SLEPC
  // Initialize SLEPc
  SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
#endif

  sub_systems_manager.petsc_initialized = true;

  // Determine if PETSc initialised MPI and is then responsible for MPI finalization
  if (!mpi_init_status && MPIinitialized())
  {
    sub_systems_manager.petsc_controls_mpi = true;
  }

#else
  error("DOLFIN has not been configured for PETSc.");
#endif

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::initZoltan()
{
#ifdef HAVE_ZOLTAN
  if ( sub_systems_manager.zoltan_initialized )
  {
    return false;
  }

  message(1, "Initializing Zoltan (ignoring command-line arguments).");

  int argc = 0;
  char** argv = NULL;
  float version;

  Zoltan_Initialize(argc, argv, &version);

  sub_systems_manager.zoltan_initialized = true;

#else
  error("DOLFIN has not been configured with Zoltan.");
#endif

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::initZoltan(int argc, char* argv[])
{
#ifdef HAVE_ZOLTAN
  if ( sub_systems_manager.zoltan_initialized )
  {
    return false;
  }

  message(1, "Initializing Zoltan with given command-line arguments.");

  float version;

  // Initialize Zoltan
  Zoltan_Initialize(argc, argv, &version);

  sub_systems_manager.petsc_initialized = true;

#else
  error("DOLFIN has not been configured with Zoltan .");
#endif

  return true;
}
//-----------------------------------------------------------------------------
void SubSystemsManager::finalizeMPI()
{
#ifdef HAVE_MPI
  if (sub_systems_manager.petsc_controls_mpi)
  {
    return;
  }

  // Finalise MPI if required
  if (MPIinitialized() && !(--mpi_init_sema_))
  {
    MPI_Finalize();
  }
#else
  // Do nothing
#endif
}
//-----------------------------------------------------------------------------
void SubSystemsManager::finalizePETSc()
{
  if (!sub_systems_manager.petsc_controls_mpi && mpi_init_sema_ > 1)
  {
    return;
  }

#ifdef HAVE_PETSC
  if (sub_systems_manager.petsc_initialized)
  {
    PetscFinalize();

#ifdef HAVE_SLEPC
    SlepcFinalize();
#endif
  }
#else
  // Do nothing
#endif
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::MPIinitialized()
{
  // This function not affected if MPI_Finalize has been called. It returns
  // true if MPI_Init has been called at any point, even if MPI_Finalize has
  // been called.

#ifdef HAVE_MPI
  int initialized;
  MPI_Initialized(&initialized);

  if (initialized) return true;
  else if (!initialized) return false;
  else
  {
    error("MPI_Initialized has returned an unknown initialization status");
    return false;
  }
#else
  // DOLFIN is not configured for MPI (it might be through PETSc)
  return false;
#endif
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
