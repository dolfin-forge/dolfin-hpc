// Copyright (C) 2007 Magnus Vikstrøm.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2008.
// Modified by Niclas Jansson, 2009-2010.
//
// First added:  2007-11-30
// Last changed: 2010-06-08

#ifndef __MPI_helper_H
#define __MPI_helper_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif 

#include <dolfin/common/types.h>

namespace dolfin
{
  /// This class provides utility functions for easy access of the number of 
  /// processes and current process number.

  class Mesh;

  class MPI
  {
  public:

    /// Return proccess number
    static uint processNumber();

    /// Return number of processes
    static uint numProcesses();

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

    /// Reorder MPI communicator
    static void reorderComm(Mesh& mesh);

#ifdef HAVE_MPI
    static MPI_Comm DOLFIN_COMM;
#else
    static int DOLFIN_COMM;
#endif

  private:
    static dolfin::real start_time;
    static int this_process, num_processes;
    static bool _this_process, _num_processes;
    static bool _dolfin_comm;
  };
}

#endif
