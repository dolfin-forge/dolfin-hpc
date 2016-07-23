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

  /// Return process number
  static uint processNumber();

  /// Return number of processes
  static uint numProcesses();

  /// Return group number
  static uint groupNumber();

  /// Return number of groups
  static uint numGroups();

  /// Return process number in world
  static uint processGlobalNumber();

  /// Return number of processes in world
  static uint numGlobalProcesses();

  /// Return if the given rank is valid
  static bool processIsValid(uint rank);

  /// Return seed value for current rank
  static uint processRandomSeed();

  /// Get process offset given number of local elements
  static void processOffset(uint local, uint& offset);

  ///
  static void numGlobalSum(uint local, uint& global);

  ///
  static void numGlobalMin(uint local, uint& global);

  ///
  static void numGlobalMax(uint local, uint& global);

  ///
  static void AllReduceSum(real local, real& global);

  ///
  static void AllReduceMin(real local, real& global);

  ///
  static void AllReduceMax(real local, real& global);

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

  /// Reorder MPI communicator
  static void reorderComm(Mesh& mesh);

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
