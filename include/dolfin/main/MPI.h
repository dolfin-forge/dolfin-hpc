// Copyright (C) 2007 Magnus Vikstrøm.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2008.
// Modified by Niclas Jansson, 2009-2010.
// Modified by Aurelien Larcher, 2016-2017.
//
// First added:  2007-11-30
// Last changed: 2017-08-24

#ifndef __DOLFIN_MPI_H
#define __DOLFIN_MPI_H

#include <dolfin/common/types.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_MPI
#define DOLFIN_COMM_NULL  MPI_COMM_NULL
#else
#define DOLFIN_COMM_NULL  0
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

  /// Return process rank in local communicator
  static uint rank();

  /// Return local communicator size
  static uint size();

  /// Return if the given rank is valid
  static bool is_valid_rank(uint rank);

  /// Return if the current process is the root process in local communicator
  static bool is_root();

  /// Return group identifier
  static uint group_id();

  /*
   *  Global communicator
   */

  /// Return process rank in global communicator
  static uint global_rank();

  /// Return global communicator size
  static uint global_size();

  /// Return if the current process is the root process in global communicator
  static bool is_global_root();

  /// Return number of groups in global communicator
  static uint num_groups();

  /*
   *  Global communicator
   */

  /// Return seed value for current process
  static uint seed();

  ///
  static
  void offset(uint xl, uint& offset, Communicator& comm = MPI::DOLFIN_COMM);

  ////
  enum ReductionType { sum, min, max };

  //// Wrap in a template function to allow use of functors
  template<class T>
  static void bcast(T* x, int n, int r, Communicator& comm = MPI::DOLFIN_COMM);

  //// Wrap in a template function to allow use of functors
  template<int R, class T>
  static void all_reduce(T x, T& r, Communicator& comm = MPI::DOLFIN_COMM);

  //// Wrap in a template function to allow use of functors
  template<class T>
  static int sendrecv(T* s, int ns, int src, T* r, int nr, int dst, int tg,
                      Communicator& comm = MPI::DOLFIN_COMM);

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
  static void finiComm();

  static Communicator DOLFIN_COMM_WORLD;
  static Communicator DOLFIN_COMM_SELF;
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

// Defines default type to MPI communication domain.

typedef MPI::Communicator Comm;

#define DOLFIN_COMM_WORLD MPI::DOLFIN_COMM_WORLD
#define DOLFIN_COMM_SELF  MPI::DOLFIN_COMM_SELF
#define DOLFIN_COMM       MPI::DOLFIN_COMM

//-----------------------------------------------------------------------------

#if HAVE_MPI

//-----------------------------------------------------------------------------
template<>
inline void MPI::bcast(uint* x, int n, int r, Communicator& comm)
{
  MPI_Bcast(x, n, MPI_UNSIGNED, r, comm);
}
//-----------------------------------------------------------------------------
template<>
inline void MPI::bcast(real* x, int n, int r, Communicator& comm)
{
  MPI_Bcast(x, n, MPI_DOUBLE, r, comm);
}
//-----------------------------------------------------------------------------
template<>
inline void MPI::all_reduce<MPI::sum>(uint x, uint& r, Communicator& comm)
{
  MPI_Allreduce(&x, &r, 1, MPI_UNSIGNED, MPI_SUM, comm);
}
//-----------------------------------------------------------------------------
template<>
inline void MPI::all_reduce<MPI::min>(uint x, uint& r, Communicator& comm)
{
  MPI_Allreduce(&x, &r, 1, MPI_UNSIGNED, MPI_MIN, comm);
}
//-----------------------------------------------------------------------------
template<>
inline void MPI::all_reduce<MPI::max>(uint x, uint& r, Communicator& comm)
{
  MPI_Allreduce(&x, &r, 1, MPI_UNSIGNED, MPI_MAX, comm);
}
//-----------------------------------------------------------------------------
template<>
inline void MPI::all_reduce<MPI::sum>(real x, real& r, Communicator& comm)
{
  MPI_Allreduce(&x, &r, 1, MPI_DOUBLE, MPI_SUM, comm);
}
//-----------------------------------------------------------------------------
template<>
inline void MPI::all_reduce<MPI::min>(real x, real& r, Communicator& comm)
{
  MPI_Allreduce(&x, &r, 1, MPI_DOUBLE, MPI_MIN, comm);
}
//-----------------------------------------------------------------------------
template<>
inline void MPI::all_reduce<MPI::max>(real x, real& r, Communicator& comm)
{
  MPI_Allreduce(&x, &r, 1, MPI_DOUBLE, MPI_MAX, comm);
}
//-----------------------------------------------------------------------------
template<>
inline int MPI::sendrecv(bool* s, int ns, int src, bool* r, int nr, int dst,
                         int tg, Communicator& comm)
{
  MPI_Status status;
  int recv_count;
  int bs = ns * sizeof(bool);
  int br = nr * sizeof(bool);
  MPI_Sendrecv(s, bs, MPI_BYTE, src, tg, r, br, MPI_BYTE, dst, tg, comm,
               &status);
  MPI_Get_count(&status, MPI_BYTE, &recv_count);
  return recv_count / sizeof(bool);
}
//-----------------------------------------------------------------------------
template<>
inline int MPI::sendrecv(int* s, int ns, int src, int* r, int nr, int dst,
                         int tg, Communicator& comm)
{
  MPI_Status status;
  int recv_count;
  MPI_Sendrecv(s, ns, MPI_INT, src, tg, r, nr, MPI_INT, dst, tg, comm, &status);
  MPI_Get_count(&status, MPI_INT, &recv_count);
  return recv_count;
}
//-----------------------------------------------------------------------------
template<>
inline int MPI::sendrecv(uint* s, int ns, int src, uint* r, int nr, int dst,
                         int tg, Communicator& comm)
{
  MPI_Status status;
  int recv_count;
  MPI_Sendrecv(s, ns, MPI_UNSIGNED, src, tg, r, nr, MPI_UNSIGNED, dst, tg, comm,
               &status);
  MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
  return recv_count;
}
//-----------------------------------------------------------------------------
template<>
inline int MPI::sendrecv(real* s, int ns, int src, real* r, int nr, int dst,
                         int tg, Communicator& comm)
{
  MPI_Status status;
  int recv_count;
  MPI_Sendrecv(s, ns, MPI_DOUBLE, src, tg, r, nr, MPI_DOUBLE, dst, tg, comm,
               &status);
  MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
  return recv_count;
}
//-----------------------------------------------------------------------------

#else

//-----------------------------------------------------------------------------
template<class T>
inline void MPI::bcast(T* x, int n, int r, Communicator& comm)
{
}

//-----------------------------------------------------------------------------
template<int R, class T>
inline void MPI::all_reduce(T x, T& r, Communicator& comm)
{ r = x; }

//-----------------------------------------------------------------------------
template<class T>
inline int MPI::sendrecv(T* s, int ns, int src, T* r, int nr, int dst, int tg,
                         Communicator& comm)
{
  return 0;
}

//-----------------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MPI_H */
