// Copyright (C) 2007 Magnus Vikstrøm.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2008.
// Modified by Niclas Jansson, 2009-2010.
// Modified by Aurelien Larcher, 2016-2017.
//
// First added:  2007-11-30
// Last changed: 2010-06-08

#ifndef __DOLFIN_MPI_H
#define __DOLFIN_MPI_H

#include <dolfin/common/types.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif


namespace dolfin
{

/// This class provides utility functions for easy access of the number of
/// processes and current process number.

class Mesh;

class MPI
{

public:

#ifdef HAVE_MPI
  typedef MPI_Comm  Communicator;
#else
  typedef int       Communicator;
#endif

  /*
   *  Local communicator
   */

  /// Return process rank in xl communicator
  static uint rank();

  /// Return xl communicator size
  static uint size();

  /// Return if the given rank is valid
  static bool is_valid_rank(uint rank);

  /// Return if the current process is the root process in xl communicator
  static bool is_root();

  /// Return group identifier
  static uint group_id();

  /*
   *  Global communicator
   */

  /// Return process rank in xg communicator
  static uint global_rank();

  /// Return xg communicator size
  static uint global_size();

  /// Return if the current process is the root process in xg communicator
  static bool is_global_root();

  /// Return number of groups in xg communicator
  static uint num_groups();

  /*
   *  Global communicator
   */

  /// Return seed value for current process
  static uint seed();

  ///
  static
  void offset(uint xl, uint& offset, Communicator& comm = MPI::DOLFIN_COMM);

  ///
  static
  void allReduceSum(uint xl, uint& xg, Communicator& comm = MPI::DOLFIN_COMM);

  ///
  static
  void allReduceMin(uint xl, uint& xg, Communicator& comm = MPI::DOLFIN_COMM);

  ///
  static
  void allReduceMax(uint xl, uint& xg, Communicator& comm = MPI::DOLFIN_COMM);

  ///
  static
  void allReduceSum(real xl, real& xg, Communicator& comm = MPI::DOLFIN_COMM);

  ///
  static
  void allReduceMin(real xl, real& xg, Communicator& comm = MPI::DOLFIN_COMM);

  ///
  static
  void allReduceMax(real xl, real& xg, Communicator& comm = MPI::DOLFIN_COMM);

  /// Start MPI timer
  static void startTimer();

  /// Start MPI timer with external counter;
  static void startTimer(dolfin::real& stime);

  /// Stop MPI timer
  static real stopTimer();

  /// Stop MPI timer
  static real stopTimer(dolfin::real& stime);

  /// Setup DOLFIN_COMM MPI communicator
  static void initComm(int ngroups = 0);

  static Communicator DOLFIN_COMM_WORLD;
  static Communicator DOLFIN_COMM;

private:

  static real time_;
  static bool init_;

  typedef struct {
    int   global_rank;
    int   global_size;
    int   group_cnt;
    int   group_idx;
    int   rank;
    int   size;
    int   seed;
  } Context;

  static Context ctx_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_MPI_H */
