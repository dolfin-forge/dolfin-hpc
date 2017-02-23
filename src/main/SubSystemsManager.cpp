// Copyright (C) 2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Niclas Jansson 2009-2015.
// Modified by Aurelien Larcher 2017.
//
// First added:  2008-01-07
// Last changed: 2017-02-23

#include <dolfin/main/SubSystemsManager.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

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

#define SUBSYSTEM_RETURN_IF_INITIALIZED(S) \
	if (sub_systems_manager.iset(S::flag)) { return false; }

#define SUBSYSTEM_RETURN_IF_CONSUMERS(S) \
  if (S::sema > 0) { return false; }

#define SUBSYSTEM_SET_INIT(S) \
  sub_systems_manager.init(S::flag);

#define SUBSYSTEM_SET_FINI(S) \
  sub_systems_manager.fini(S::flag);

#define SUBSYSTEM_INITIALIZED(S) \
  sub_systems_manager.iset(S::flag)

#define SUBSYSTEM_ERROR_NOT_ENABLED(S) \
  error("DOLFIN has not been configured with " #S ".")

#define SUBSYSTEM_CHECK_SEMA(S) \
	if (S::sema == 0) \
  { \
    error("SubsystemsManager : finalization but " #S " has no consumer"); \
  }

//--- STATIC ------------------------------------------------------------------

SubSystemsManager SubSystemsManager::sub_systems_manager;

int SubSystemsManager::MPI::sema    = 0;
int SubSystemsManager::PETSc::sema  = 0;
int SubSystemsManager::Zoltan::sema = 0;

//-----------------------------------------------------------------------------
SubSystemsManager::SubSystemsManager() :
    count_(0),
    state_(0)
{
}
//-----------------------------------------------------------------------------
SubSystemsManager::SubSystemsManager(SubSystemsManager const& other) :
    count_(0),
    state_(0)
{
  error("Should not be using copy constructor of SubSystemsManager.");
}
//-----------------------------------------------------------------------------
SubSystemsManager::~SubSystemsManager()
{
  SubSystemsManager::fini();
}
//-----------------------------------------------------------------------------
int SubSystemsManager::init(int argc, char* argv[], uint n)
{
  if (count_ == 0 )
  {
#ifdef HAVE_MPI
    SubSystemsManager::MPI::init(argc, argv, n);
#endif

#ifdef HAVE_PETSC
    SubSystemsManager::PETSc::init(argc, argv);
#endif

#ifdef HAVE_ZOLTAN
    SubSystemsManager::Zoltan::init(argc, argv);
#endif
  }
  return ++count_;
}
//-----------------------------------------------------------------------------
int SubSystemsManager::fini()
{
  if (count_ > 0 )
  {
    // Finalize subsystems in the correct order
#ifdef HAVE_ZOLTAN
    SubSystemsManager::Zoltan::fini();
#endif

#ifdef HAVE_PETSC
    SubSystemsManager::PETSc::fini();
#endif

#ifdef HAVE_MPI
    SubSystemsManager::MPI::fini();
#endif

    count_ = 0;
  }
  return count_;
}
//-------------------------------------------------------------------------
void SubSystemsManager::disp() const
{
  section("SubSystemsManager");
  message("count : %u", count_);
  message("state : %u", state_);
  message("      - MPI      %u", iset(mpi));
  message("      - PETSc    %u", iset(petsc));
  message("      - PETScMPI %u", iset(petscmpi));
  message("      - JANPACK  %u", iset(janpack));
  message("      - Zoltan   %u", iset(zoltan));
  message("      - SLEPc    %u", iset(slepc));
  endblock();
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::MPI::init(int argc, char* argv[], uint n)
{
#ifdef HAVE_MPI

  ++MPI::sema;
  SUBSYSTEM_RETURN_IF_INITIALIZED(MPI);

  /// MPI was initialized by another subsystem: sema is incremented to avoid
  /// finalization when all consumers in DOLFIN have returned.
  if (MPI::initialized()) { ++MPI::sema; return false; }

#ifdef HAVE_JANPACK_MPI
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  SUBSYSTEM_SET_INIT(JANPACK);
#else
  MPI_Init(&argc, &argv);
#endif /* HAVE_JANPACK_MPI */

  dolfin::MPI::initComm(n);
  SUBSYSTEM_SET_INIT(MPI);

#else

  SUBSYSTEM_ERROR_NOT_ENABLED(MPI);

#endif /* HAVE_MPI */

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::MPI::fini()
{
#ifdef HAVE_MPI

  SUBSYSTEM_CHECK_SEMA(MPI);
  --MPI::sema;
  SUBSYSTEM_RETURN_IF_CONSUMERS(MPI)

  /// PETSc is responsible for PETSc finalize
  SUBSYSTEM_RETURN_IF_INITIALIZED(PETScMPI);

  if (!MPI::initialized())
  {
    error("SubsystemsManager : finalization called on uninitialized MPI");
  }

  dolfin::MPI::finiComm();
  MPI_Finalize();
  SUBSYSTEM_SET_FINI(MPI);

#else

  SUBSYSTEM_ERROR_NOT_ENABLED(MPI);

#endif /* HAVE_MPI */
  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::MPI::initialized()
{
#ifdef HAVE_MPI

  int initialized;
  MPI_Initialized(&initialized);
  return (initialized > 0);

#else

  SUBSYSTEM_ERROR_NOT_ENABLED(MPI);
  return false;

#endif /* HAVE_MPI */
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::PETSc::init(int argc, char* argv[])
{
#ifdef HAVE_PETSC

  ++PETSc::sema;
  SUBSYSTEM_RETURN_IF_INITIALIZED(PETSc);

  // Get status of MPI before PETSc initialization
  bool const mpi_init_status = MPI::initialized();

  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  SUBSYSTEM_SET_INIT(PETSc);

  // If PETSc initialized MPI, it is responsible for MPI finalization
  if (!mpi_init_status && MPI::initialized())
  {
    SUBSYSTEM_SET_INIT(PETScMPI);
  }

#ifdef HAVE_SLEPC

  // Initialize SLEPc
  SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  SUBSYSTEM_SET_INIT(SLEPc);

#endif /*  HAVE_SLEPC */

#else

  SUBSYSTEM_ERROR_NOT_ENABLED(PETSc);

#endif /* HAVE_PETSC */

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::PETSc::fini()
{
#ifdef HAVE_PETSC

  SUBSYSTEM_CHECK_SEMA(PETSc);
  --PETSc::sema;
  SUBSYSTEM_RETURN_IF_CONSUMERS(PETSc)

  /// PETSc is responsible for MPI and there are still consumers
  if (SUBSYSTEM_INITIALIZED(PETScMPI) && (MPI::sema > 1)) { return false; }

#ifdef HAVE_SLEPC
    SlepcFinalize();
    SUBSYSTEM_SET_FINI(SLEPc);
#endif

  PetscFinalize();
  if (SUBSYSTEM_INITIALIZED(PETScMPI)) SUBSYSTEM_SET_FINI(PETScMPI);
  SUBSYSTEM_SET_FINI(PETSc);

#else

  SUBSYSTEM_ERROR_NOT_ENABLED(PETSc);

#endif /* HAVE_PETSC */

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::Zoltan::init(int argc, char* argv[])
{
#ifdef HAVE_ZOLTAN

  ++Zoltan::sema;
  SUBSYSTEM_RETURN_IF_INITIALIZED(Zoltan);

  // Initialize Zoltan
  float version;
  Zoltan_Initialize(argc, argv, &version);

  SUBSYSTEM_SET_INIT(Zoltan)

#else

  SUBSYSTEM_ERROR_NOT_ENABLED(Zoltan);

#endif /* HAVE_ZOLTAN */

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::Zoltan::fini()
{
#ifdef HAVE_ZOLTAN

  SUBSYSTEM_CHECK_SEMA(Zoltan);
  --Zoltan::sema;
  SUBSYSTEM_RETURN_IF_CONSUMERS(Zoltan)
  SUBSYSTEM_SET_FINI(Zoltan);

#else

  SUBSYSTEM_ERROR_NOT_ENABLED(Zoltan);

#endif /* HAVE_ZOLTAN */

  return true;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
