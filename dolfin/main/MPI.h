// Copyright (C) 2007 Magnus Vikstrøm.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2008.
//
// First added:  2007-11-30
// Last changed: 2008-01-07

#ifndef __MPI_helper_H
#define __MPI_helper_H

#ifdef HAS_MPI
#include <mpi.h>
#endif 

#include <dolfin/common/types.h>

namespace dolfin
{
  /// This class provides utility functions for easy access of the number of 
  /// processes and current process number.
  
  class MPI
  {
  public:

    /// Return proccess number
    static uint processNumber();

    /// Return number of processes
    static uint numProcesses();

    /// Determine whether we should broadcast (based on current parallel policy)
    static bool broadcast();

    /// Determine whether we should receive (based on current parallel policy)
    static bool receive();

    /// Start MPI timer
    static void startTimer();
    
    /// Start MPI timer with external counter;
    static void startTimer(dolfin::real& stime);

    /// Stop MPI timer
    static real stopTimer();

    /// Stop MPI timer
    static real stopTimer(dolfin::real& stime);

    /// Setup DOLFIN_COMM MPI communicator
    static void initComm();

    //#ifdef HAS_MPI
    static MPI_Comm DOLFIN_COMM;
    //#endif

  private:
    static dolfin::real start_time;
    static int this_process, num_processes;
    static bool _this_process, _num_processes;
    static bool _dolfin_comm;
  };
}

#endif
