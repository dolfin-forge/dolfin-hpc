// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/main/init.h>

#include <dolfin/config/dolfin_config.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/main/SubSystemsManager.h>

#include <cstdlib>

namespace dolfin
{

//-----------------------------------------------------------------------------
void dolfin_init( int    argc,
                  char * argv[],
                  long   wallclock_limit,
                  int    parallel_groups )
{

  //--- Initialize subsystems and print banner
  char const * verbosity = std::getenv( "DOLFIN_VERBOSE" );
  if ( verbosity != NULL )
  {
    verbose( std::atoi( verbosity ) );
    message( 1, "Set verbosity level to %i", std::atoi( verbosity ) );
  }

  //--- Initialize subsystems and print banner
  if ( SubSystemsManager::start( argc, argv, parallel_groups, wallclock_limit )
       == 1 )
  {
#ifdef DOLFIN_HAVE_MPI
    if ( MPI::global_rank() == 0 )
    {
      message( "Initializing DOLFIN version %s\n%s\n\nRunning on %d %s (%u %s)",
               DOLFIN_VERSION, DOLFIN_BUILD_INFO,
               MPI::global_size(),
               MPI::global_size() > 1 ? "processes" : "process",
               MPI::num_groups(),
               MPI::num_groups() > 1 ? "groups" : "group" );
    }
    else
    {
      silence();
    }
#else
    message( "Initializing DOLFIN version %s\n%s",
             DOLFIN_VERSION, DOLFIN_BUILD_INFO );
#endif
  }
}
//-----------------------------------------------------------------------------
void dolfin_finalize()
{
}
//-----------------------------------------------------------------------------

} // end namespace dolfin
