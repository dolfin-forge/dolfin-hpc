// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/main/init.h>

#include <dolfin/common/constants.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/main/SubSystemsManager.h>

#include <cstdlib>
#include <cstring>
#include <unistd.h>

//-----------------------------------------------------------------------------
void dolfin::dolfin_init(int argc, char * argv[], long w_limit, int n)
{
  //--- Initialize subsystems and print banner

  int init_count = SubSystemsManager::start(argc, argv, n, w_limit);
  if (init_count == 1)
  {
#ifdef HAVE_MPI
    if (MPI::global_rank() == 0)
    {
      message("Initializing DOLFIN version %s\n%s\n\nRunning on %d %s (%u %s)",
              DOLFIN_VERSION,
              DOLFIN_BUILD_INFO,
              dolfin::MPI::global_size(),
              dolfin::MPI::global_size() > 1 ? "processes" : "process",
              dolfin::MPI::num_groups(),
              dolfin::MPI::num_groups()  > 1 ? "groups"    : "group");
    }
    else
    {
      silence();
    }
#else
    message("Initializing DOLFIN version %s\n%s",
            DOLFIN_VERSION, DOLFIN_BUILD_INFO);
#endif
  }

}
//-----------------------------------------------------------------------------
void dolfin::dolfin_finalize()
{
}
//-----------------------------------------------------------------------------
