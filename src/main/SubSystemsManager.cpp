// Copyright (C) 2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/main/SubSystemsManager.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/maybe_unused.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

#ifdef HAVE_PETSC
#include <petsc.h>
#endif

#ifdef HAVE_ZOLTAN
#include <zoltan_cpp.h>
#endif

#ifdef HAVE_TRILINOS
#include <Kokkos_Core.hpp>
#include <Epetra_config.h>

#ifdef DOLFIN_HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#endif

#include <cstdlib>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------

int SubSystemsManager::MPI::sema      = 0;
int SubSystemsManager::PETSc::sema    = 0;
int SubSystemsManager::Zoltan::sema   = 0;
int SubSystemsManager::Trilinos::sema = 0;

//-----------------------------------------------------------------------------
SubSystemsManager::SubSystemsManager()
  : count_( 0 )
  , state_( 0 )
{
}
//-----------------------------------------------------------------------------
SubSystemsManager::SubSystemsManager( SubSystemsManager const & )
  : count_( 0 )
  , state_( 0 )
{
  error( "Should not be using copy constructor of SubSystemsManager." );
}
//-----------------------------------------------------------------------------
SubSystemsManager::~SubSystemsManager()
{
  SubSystemsManager::finalize();
}
//-----------------------------------------------------------------------------
int SubSystemsManager::initialize( int    argc,
                                   char * argv[],
                                   uint   n,
                                   long   w_limit )
{
  if ( count_ == 0 )
  {
    char const * verbosity = std::getenv( "DOLFIN_VERBOSE" );
    if ( verbosity != nullptr )
    {
      verbose( std::atoi( verbosity ) );
    }

#ifdef DOLFIN_HAVE_MPI
    SubSystemsManager::MPI::initialize( argc, argv, n );
#else
    MAYBE_UNUSED( n );
#endif

#ifdef HAVE_PETSC
    SubSystemsManager::PETSc::initialize( argc, argv );
#endif

#ifdef HAVE_ZOLTAN
    SubSystemsManager::Zoltan::initialize( argc, argv );
#endif

#if !defined( DOLFIN_HAVE_MPI ) and !defined( HAVE_PETSC ) and !defined( HAVE_ZOLTAN )
    MAYBE_UNUSED( argc );
    MAYBE_UNUSED( argv );
#endif

#ifdef HAVE_TRILINOS
    SubSystemsManager::Trilinos::initialize( argc, argv );
#endif

    // Set wall clock limit
    timer_.set_limit( w_limit );
  }

  if ( verbose() > 0 )
    this->disp();

  return ++count_;
}
//-----------------------------------------------------------------------------
int SubSystemsManager::finalize()
{
  if ( count_ > 0 )
  {
    // Finalize subsystems in the correct order
#ifdef HAVE_ZOLTAN
    SubSystemsManager::Zoltan::finalize();
#endif

#ifdef HAVE_PETSC
    SubSystemsManager::PETSc::finalize();
#endif

#ifdef HAVE_TRILINOS
    SubSystemsManager::Trilinos::finalize();
#endif

#ifdef DOLFIN_HAVE_MPI
    SubSystemsManager::MPI::finalize();
#endif

    count_ = 0;
  }
  return count_;
}
//-------------------------------------------------------------------------
void SubSystemsManager::disp() const
{
  section( "SubSystemsManager" );
  message( "count : %u", count_ );
  message( "state : %u", state_ );
  message( "      - MPI      %u", iset( mpi ) );
  message( "      - PETSc    %u", iset( petsc ) );
  message( "      - PETScMPI %u", iset( petscmpi ) );
  message( "      - JANPACK  %u", iset( janpack ) );
  message( "      - Zoltan   %u", iset( zoltan ) );
  message( "      - Trilinos %u", iset( trilinos ) );
  end();
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::MPI::initialize( int argc, char * argv[], uint n )
{
#ifdef DOLFIN_HAVE_MPI

  ++MPI::sema;
  if ( SubSystemsManager::instance().iset( MPI::flag ) )
  {
    return false;
  }

  /// MPI was initialized by another subsystem: sema is incremented to avoid
  /// finalization when all consumers in DOLFIN have returned.
  if ( MPI::initialized() )
  {
    ++MPI::sema;
    return false;
  }

#ifdef HAVE_JANPACK_MPI
  int provided;
  dolfin::MPI::check_error(
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &provided ) );
  SubSystemsManager::instance().initialize( JANPACK::flag );
#else
  dolfin::MPI::check_error( MPI_Init( &argc, &argv ) );
#endif /* HAVE_JANPACK_MPI */

  dolfin::MPI::initComm( n );
  SubSystemsManager::instance().initialize( MPI::flag );

#else

  error( "DOLFIN has not been configured with MPI." );

  MAYBE_UNUSED( argc )
  MAYBE_UNUSED( argv )
  MAYBE_UNUSED( n )
#endif /* DOLFIN_HAVE_MPI */

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::MPI::finalize()
{
#ifdef DOLFIN_HAVE_MPI
  if ( MPI::sema == 0 )
  {
    error( "SubsystemsManager : finalization but MPI has no consumer" );
  }

  --MPI::sema;
  if ( MPI::sema > 0 )
  {
    return false;
  }

  /// PETSc is responsible for PETSc finalize
  if ( SubSystemsManager::instance().iset( PETScMPI::flag ) )
  {
    return false;
  }

  if ( !MPI::initialized() )
  {
    error( "SubsystemsManager : finalization called on uninitialized MPI" );
  }

  dolfin::MPI::finiComm();
  dolfin::MPI::check_error( MPI_Finalize() );
  SubSystemsManager::instance().finalize( MPI::flag );

#else

  error( "DOLFIN has not been configured with MPI." );

#endif /* DOLFIN_HAVE_MPI */
  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::MPI::initialized()
{
#ifdef DOLFIN_HAVE_MPI

  int initialized;
  dolfin::MPI::check_error( MPI_Initialized( &initialized ) );
  return ( initialized > 0 );

#else

  error( "DOLFIN has not been configured with MPI." ) return false;

#endif /* DOLFIN_HAVE_MPI */
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::PETSc::initialize( int argc, char * argv[] )
{
#ifdef HAVE_PETSC

  ++PETSc::sema;

  if ( SubSystemsManager::instance().iset( PETSc::flag ) )
  {
    return false;
  }

#ifdef DOLFIN_HAVE_MPI
  // Get status of MPI before PETSc initialization
  bool const mpi_init_status = MPI::initialized();
#endif

  // Initialize PETSc
  PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
  SubSystemsManager::instance().initialize( PETSc::flag );

  // remove PETSc Signal handling
  PetscPopSignalHandler();

#ifdef DOLFIN_HAVE_MPI
  // If PETSc initialized MPI, it is responsible for MPI finalization
  if ( !mpi_init_status and MPI::initialized() )
  {
    SubSystemsManager::instance().initialize( PETScMPI::flag );
  }
#endif

#else

  error( "DOLFIN has not been configured with PETSc." );

#endif /* HAVE_PETSC */

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::PETSc::finalize()
{
#ifdef HAVE_PETSC

  if ( PETSc::sema == 0 )
  {
    error( "SubsystemsManager : finalization but PETSc has no consumer" );
  }

  --PETSc::sema;

  if ( PETSc::sema > 0 )
  {
    return false;
  }

#ifdef DOLFIN_HAVE_MPI
  /// PETSc is responsible for MPI and there are still consumers
  if ( SubSystemsManager::instance().iset( PETScMPI::flag )
       and ( MPI::sema > 1 ) )
  {
    return false;
  }
#endif

  PetscFinalize();

#ifdef DOLFIN_HAVE_MPI
  if ( SubSystemsManager::instance().iset( PETScMPI::flag ) )
    SubSystemsManager::instance().finalize( PETScMPI::flag );
#endif

  SubSystemsManager::instance().finalize( PETSc::flag );

#else

  error( "DOLFIN has not been configured with PETSc." );

#endif /* HAVE_PETSC */

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::Zoltan::initialize( int argc, char * argv[] )
{
#ifdef HAVE_ZOLTAN

  ++Zoltan::sema;

  if ( SubSystemsManager::instance().iset( Zoltan::flag ) )
  {
    return false;
  }

  // Initialize Zoltan
  float version;
  Zoltan_Initialize( argc, argv, &version );

  SubSystemsManager::instance().initialize( Zoltan::flag );

#else

  error( "DOLFIN has not been configured with Zoltan." );
  MAYBE_UNUSED( argc )
  MAYBE_UNUSED( argv )

#endif /* HAVE_ZOLTAN */

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::Zoltan::finalize()
{
#ifdef HAVE_ZOLTAN

  if ( Zoltan::sema == 0 )
  {
    error( "SubsystemsManager : finalization but Zoltan has no consumer" );
  }

  --Zoltan::sema;

  if ( Zoltan::sema > 0 )
  {
    return false;
  }

  SubSystemsManager::instance().finalize( Zoltan::flag );

#else

  error( "DOLFIN has not been configured with Zoltan." );

#endif /* HAVE_ZOLTAN */

  return true;
} //-----------------------------------------------------------------------------
bool SubSystemsManager::Trilinos::initialize( int argc, char * argv[] )
{
#ifdef HAVE_TRILINOS
  ++Trilinos::sema;

  if ( SubSystemsManager::instance().iset( Trilinos::flag ) )
  {
    return false;
  }

#ifdef DOLFIN_HAVE_MPI
  // Trilinos has to be initialized after MPI
  if ( not MPI::initialized() )
    error( "SubsystemsManager : Trilinos has to be initialized after MPI" );
  else
    SubSystemsManager::instance().initialize( TrilinosMPI::flag );
#endif

  // Initialize Trilinos
  Kokkos::initialize( argc, argv );
  SubSystemsManager::instance().initialize( Trilinos::flag );

#else

  error( "DOLFIN has not been configured with Trilinos." );

#endif /* HAVE_TRILINOS */

  return true;
}
//-----------------------------------------------------------------------------
bool SubSystemsManager::Trilinos::finalize()
{
#ifdef HAVE_TRILINOS

  if ( Trilinos::sema == 0 )
  {
    error( "SubsystemsManager : finalization but Trilinos has no consumer" );
  }

  --Trilinos::sema;

  if ( Trilinos::sema > 0 )
  {
    return false;
  }

#ifdef DOLFIN_HAVE_MPI
  if ( SubSystemsManager::instance().iset( TrilinosMPI::flag )
       and ( MPI::sema > 1 ) )
  {
    error( "SubsystemsManager : Trilinos has to be finalized before MPI" );
  }
#endif

  // finalize trilinos
  Kokkos::finalize();
  SubSystemsManager::instance().finalize( Trilinos::flag );

#else

  error( "DOLFIN has not been configured with Trilinos." );

#endif /* HAVE_TRILINOS */

  return true;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
