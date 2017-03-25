// Copyright (C) 2007 Magnus Vikstrøm.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007, 2008.
// Modified by Anders Logg, 2007.
// Modified by Niclas Jansson, 2009-2010.
// Modified by Aurelien Larcher, 2015-2017.
//
// First added:  2007-11-30
// Last changed: 2017-02-17

#include <dolfin/main/MPI.h>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/SubSystemsManager.h>

#include <cstring>
#include <ctime>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------
MPI::Communicator MPI::DOLFIN_COMM_WORLD;
MPI::Communicator MPI::DOLFIN_COMM_SELF;
MPI::Communicator MPI::DOLFIN_COMM;

real MPI::time_ = 0.0;
bool MPI::init_ = false;
MPI::Context MPI::ctx_;

#if HAVE_MPI
#define DOLFIN_MPI_SUBSYSTEM_INIT \
  if (!init_) { SubSystemsManager::MPI::init(); }
#endif

#define DOLFIN_MPI_UNIMPLEMENTED \
	error("Unimplemented without MPI support");

//-----------------------------------------------------------------------------
uint MPI::rank()
{
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.rank);
}
//-----------------------------------------------------------------------------
uint MPI::size()
{
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.size);
}
//-----------------------------------------------------------------------------
uint MPI::group_id()
{
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.group_idx);
}
//-----------------------------------------------------------------------------
uint MPI::num_groups()
{
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.group_cnt);
}
//-----------------------------------------------------------------------------
uint MPI::global_rank()
{
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.global_rank);
}
//-----------------------------------------------------------------------------
uint MPI::global_size()
{
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.global_size);
}
//-----------------------------------------------------------------------------
void MPI::startTimer()
{
#ifdef HAVE_MPI
  MPI_Barrier(MPI::DOLFIN_COMM);
  time_ = MPI_Wtime();
#else
  DOLFIN_MPI_UNIMPLEMENTED
#endif
}
//-----------------------------------------------------------------------------
real MPI::stopTimer()
{
#ifdef HAVE_MPI
  MPI_Barrier(MPI::DOLFIN_COMM);
  return (MPI_Wtime() - time_);
#else
  DOLFIN_MPI_UNIMPLEMENTED
#endif
}
//-----------------------------------------------------------------------------
void MPI::startTimer(real& stime)
{
#ifdef HAVE_MPI
  MPI_Barrier(MPI::DOLFIN_COMM);
  stime = MPI_Wtime();
#else
  DOLFIN_MPI_UNIMPLEMENTED
#endif
}
//-----------------------------------------------------------------------------
real MPI::stopTimer(real& stime)
{
#ifdef HAVE_MPI
  MPI_Barrier(MPI::DOLFIN_COMM);
  return (MPI_Wtime() - stime);
#else
  DOLFIN_MPI_UNIMPLEMENTED
#endif
}
//-----------------------------------------------------------------------------
void MPI::initComm(int ngroups)
{
  if (init_) { return; }

#ifdef HAVE_MPI

  // Initialize world
  MPI_Comm_dup(MPI_COMM_WORLD, &MPI::DOLFIN_COMM_WORLD);
  MPI_Comm_rank(MPI::DOLFIN_COMM_WORLD, &ctx_.global_rank);
  MPI_Comm_size(MPI::DOLFIN_COMM_WORLD, &ctx_.global_size);
  MPI_Comm_dup(MPI_COMM_SELF, &MPI::DOLFIN_COMM_SELF);
  ctx_.seed = std::time(0) + ctx_.global_rank;

  // Initialize group(s)
  int const wsize = ctx_.global_size;
  int const wrank = ctx_.global_rank;
  if ((ngroups > 1) && (wsize >= ngroups))
  {
    MPI_Group wgroup;
    MPI_Comm_group(MPI::DOLFIN_COMM_WORLD, &wgroup);

    int const p = wsize / ngroups;
    int const r = wsize % ngroups;
    int const k = wrank / (r ? p + 1 : p);
    int const ssize = (wrank < (p + 1) * r ? p + 1 : p);
    int const sroot = k * p + std::min(k, r);
    int srange[3] = { sroot, sroot + ssize - 1, 1 };

    MPI_Group sgroup;
    MPI_Group_range_incl(wgroup, 1, &srange, &sgroup);
    MPI_Comm_create(MPI::DOLFIN_COMM_WORLD, sgroup, &MPI::DOLFIN_COMM);
    MPI_Group_rank(sgroup, &ctx_.rank);
    MPI_Group_size(sgroup, &ctx_.size);
    ctx_.group_idx = k;
    ctx_.group_cnt = ngroups;
  }
  else
  {
    MPI_Comm_dup(MPI::DOLFIN_COMM_WORLD, &MPI::DOLFIN_COMM);
    MPI_Comm_rank(MPI::DOLFIN_COMM, &ctx_.rank);
    MPI_Comm_size(MPI::DOLFIN_COMM, &ctx_.size);
    ctx_.group_idx = 0;
    ctx_.group_cnt = 1;
  }
#endif

  init_ = true;
}
//-----------------------------------------------------------------------------
void MPI::finiComm()
{
  MPI_Comm_free(&MPI::DOLFIN_COMM);
  MPI_Comm_free(&MPI::DOLFIN_COMM_SELF);
  MPI_Comm_free(&MPI::DOLFIN_COMM_WORLD);
}
//-----------------------------------------------------------------------------
uint MPI::seed()
{
  return ctx_.seed;
}
//-----------------------------------------------------------------------------
bool MPI::is_valid_rank(uint rank)
{
  return (rank < static_cast<uint>(ctx_.size));
}
//-----------------------------------------------------------------------------
bool MPI::is_root()
{
  return (ctx_.rank == 0);
}
//-----------------------------------------------------------------------------
void MPI::offset(uint local, uint& offset, Communicator& comm)
{
  // Fool-proof as the value for rank 0 is undefined according to MPI specs
  offset = 0;
#if ( MPI_VERSION > 1 )
  MPI_Exscan(&local, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
  MPI_Scan(&local, &offset, 1, MPI_UNSIGNED, MPI_SUM, comm);
  offset -= local;
#endif
}
//-----------------------------------------------------------------------------
void MPI::allReduceSum(uint local, uint& global, Communicator& comm)
{
#if HAVE_MPI
    MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED, MPI_SUM, comm);
#else
    global = local;
#endif
}
//-----------------------------------------------------------------------------
void MPI::allReduceMin(uint local, uint& global, Communicator& comm)
{
#if HAVE_MPI
    MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED, MPI_MIN, comm);
#else
    global = local;
#endif
}
//-----------------------------------------------------------------------------
void MPI::allReduceMax(uint local, uint& global, Communicator& comm)
{
#if HAVE_MPI
    MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED, MPI_MAX, comm);
#else
    global = local;
#endif
}
//-----------------------------------------------------------------------------
void MPI::allReduceSum(real local, real& global, Communicator& comm)
{
#if HAVE_MPI
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, comm);
#else
    global = local;
#endif
}
//-----------------------------------------------------------------------------
void MPI::allReduceMin(real local, real& global, Communicator& comm)
{
#if HAVE_MPI
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MIN, comm);
#else
    global = local;
#endif
}
//-----------------------------------------------------------------------------
void MPI::allReduceMax(real local, real& global, Communicator& comm)
{
#if HAVE_MPI
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, comm);
#else
    global = local;
#endif
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
