// Copyright (C) 2007 Magnus Vikstrøm.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2008.
// Modified by Niclas Jansson, 2009-2010.
//
// First added:  2007-11-30
// Last changed: 2010-06-08

#ifndef __DOLFIN_MPI_H
#define __DOLFIN_MPI_H

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

  /*
   *  Local communicator
   */

  /// Return process rank in local communicator
  static uint rank();

  /// Return local communicator size
  static uint size();

  /// Return if the given rank is valid
  static bool is_valid_rank(uint rank);

  /// Return if the current process is the master process
  static bool is_root();

  /*
   *  Global communicator
   */

  /// Return process rank in global communicator
  static uint global_rank();

  /// Return global communicator size
  static uint global_size();

  /// Return group number
  static uint groupNumber();

  /// Return number of groups
  static uint numGroups();

  /*
   *  Global communicator
   */

  /// Return seed value for current process
  static uint seed();

  /// Get process offset given number of local elements
  static void offset(uint local, uint& offset, MPI_Comm& comm = MPI::DOLFIN_COMM);

  ///
  static void allReduceSum(uint local, uint& global, MPI_Comm& comm = MPI::DOLFIN_COMM);

  ///
  static void allReduceMin(uint local, uint& global, MPI_Comm& comm = MPI::DOLFIN_COMM);

  ///
  static void allReduceMax(uint local, uint& global, MPI_Comm& comm = MPI::DOLFIN_COMM);

  ///
  static void allReduceSum(real local, real& global, MPI_Comm& comm = MPI::DOLFIN_COMM);

  ///
  static void allReduceMin(real local, real& global, MPI_Comm& comm = MPI::DOLFIN_COMM);

  ///
  static void allReduceMax(real local, real& global, MPI_Comm& comm = MPI::DOLFIN_COMM);

  /// Start MPI timer
  static void startTimer();

  /// Start MPI timer with external counter;
  static void startTimer(dolfin::real& stime);

  /// Stop MPI timer
  static real stopTimer();

  /// Stop MPI timer
  static real stopTimer(dolfin::real& stime);

  /// Setup DOLFIN_COMM MPI communicator
  static void initComm(int n = 0);

#ifdef HAVE_MPI
  static MPI_Comm DOLFIN_COMM_WORLD;
  static MPI_Comm DOLFIN_COMM;
#else
  static int DOLFIN_COMM_WORLD;
  static int DOLFIN_COMM;
#endif

private:

  static dolfin::real start_time;
  static int this_process_world;
  static int num_processes_world;
  static int this_group;
  static int num_groups;
  static int this_process;
  static int num_processes;
  static int this_seed;
  static bool _dolfin_comm;
};

}

#endif
