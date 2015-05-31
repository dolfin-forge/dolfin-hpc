// Copyright (C) 2007 Magnus Vikstrøm.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007, 2008.
// Modified by Anders Logg, 2007.
// Modified by Niclas Jansson, 2009-2010.
// Modified by Aurelien Larcher, 2015.
//
// First added:  2007-11-30
// Last changed: 2010-01-13

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/SubSystemsManager.h>
#include <dolfin/main/MPI.h>
#include <cstring>
#include <ctime>

//-----------------------------------------------------------------------------
dolfin::real dolfin::MPI::start_time = 0.0;
int dolfin::MPI::this_process_world = 0;
int dolfin::MPI::num_processes_world = 0;
int dolfin::MPI::this_group = 0;
int dolfin::MPI::num_groups = 0;
int dolfin::MPI::this_process = 0;
int dolfin::MPI::num_processes = 0;
int dolfin::MPI::this_seed = 0;
bool dolfin::MPI::_dolfin_comm = false;
//-----------------------------------------------------------------------------
#ifdef HAVE_MPI
MPI_Comm dolfin::MPI::DOLFIN_COMM_WORLD;
MPI_Comm dolfin::MPI::DOLFIN_COMM;
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::processNumber()
{
  if (!_dolfin_comm)
  {
    SubSystemsManager::initMPI();
  }

  return static_cast<uint>(this_process);
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::numProcesses()
{
  if (!_dolfin_comm)
  {
    SubSystemsManager::initMPI();
  }

  return static_cast<uint>(num_processes);
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::groupNumber()
{
  if (!_dolfin_comm)
  {
    SubSystemsManager::initMPI();
  }

  return static_cast<uint>(this_group);
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::numGroups()
{
  if (!_dolfin_comm)
  {
    SubSystemsManager::initMPI();
  }

  return static_cast<uint>(num_groups);
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::processGlobalNumber()
{
  if (!_dolfin_comm)
  {
    SubSystemsManager::initMPI();
  }

  return static_cast<uint>(this_process_world);
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::numGlobalProcesses()
{
  if (!_dolfin_comm)
  {
    SubSystemsManager::initMPI();
  }

  return static_cast<uint>(num_processes_world);
}
#else

//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::processNumber()
{
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::numProcesses()
{
  return 1;
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::groupNumber()
{
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::numGroups()
{
  return 1;
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::processGlobalNumber()
{
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::numGlobalProcesses()
{
  return 1;
}
//-----------------------------------------------------------------------------

#endif
//-----------------------------------------------------------------------------
void dolfin::MPI::startTimer()
{
#ifdef HAVE_MPI
  MPI_Barrier(MPI::DOLFIN_COMM);
  start_time = MPI_Wtime();
#else
  error("Unimplemented without MPI support");
#endif
}
//-----------------------------------------------------------------------------
dolfin::real dolfin::MPI::stopTimer()
{
#ifdef HAVE_MPI
  MPI_Barrier(MPI::DOLFIN_COMM);
  return (MPI_Wtime() - start_time);
#else
  error("Unimplemented without MPI support");
#endif
}
//-----------------------------------------------------------------------------
void dolfin::MPI::startTimer(real& stime)
{
#ifdef HAVE_MPI
  MPI_Barrier(MPI::DOLFIN_COMM);
  stime = MPI_Wtime();
#else
  error("Unimplemented without MPI support");
#endif
}
//-----------------------------------------------------------------------------
dolfin::real dolfin::MPI::stopTimer(real& stime)
{
#ifdef HAVE_MPI
  MPI_Barrier(MPI::DOLFIN_COMM);
  return (MPI_Wtime() - stime);
#else
  error("Unimplemented without MPI support");
#endif
}
//-----------------------------------------------------------------------------
void dolfin::MPI::initComm(int n)
{
  if (_dolfin_comm)
  {
    return;
  }
#ifdef HAVE_MPI
  // Initialize world
  MPI_Comm_dup(MPI_COMM_WORLD, &MPI::DOLFIN_COMM_WORLD);
  MPI_Comm_rank(MPI::DOLFIN_COMM_WORLD, &this_process_world);
  MPI_Comm_size(MPI::DOLFIN_COMM_WORLD, &num_processes_world);
  this_seed = std::time(0) + this_process_world;

  // Initialize group(s)
  int glob_numprocs;
  MPI_Comm_size(MPI::DOLFIN_COMM_WORLD, &glob_numprocs);
  if ((n > 1) && (glob_numprocs >= n))
  {
    MPI_Group glb_group;
    MPI_Comm_group(MPI::DOLFIN_COMM_WORLD, &glb_group);
    int rank;
    MPI_Comm_rank(MPI::DOLFIN_COMM_WORLD, &rank);

    int const p = glob_numprocs / n;
    int const k = rank / (glob_numprocs % n ? p + 1 : p);
    int numprocs = (rank < (p + 1) * (glob_numprocs % n) ? p + 1 : p);
    int offset = k * p + std::min(k, glob_numprocs % n);

    // Check consistency
    //int totalprocs;
    //int value = (rank == offset ? numprocs : 0);
    //MPI_Allreduce(&value, &totalprocs, 1, MPI_INT, MPI_SUM, MPI::DOLFIN_COMM_WORLD );

    MPI_Group sub_group;
    int range[3] = { offset, offset + numprocs - 1, 1 };
    MPI_Group_range_incl(glb_group, 1, &range, &sub_group);
    MPI_Comm_create(MPI::DOLFIN_COMM_WORLD, sub_group, &MPI::DOLFIN_COMM);
    MPI_Group_rank(sub_group, &this_process);
    MPI_Group_size(sub_group, &num_processes);
    this_group = k;
    num_groups = n;
  }
  else
  {
    MPI_Comm_dup(MPI::DOLFIN_COMM_WORLD, &MPI::DOLFIN_COMM);
    MPI_Comm_rank(MPI::DOLFIN_COMM, &this_process);
    MPI_Comm_size(MPI::DOLFIN_COMM, &num_processes);
    this_group = 0;
    num_groups = 1;
  }
#endif
  _dolfin_comm = true;
}
//-----------------------------------------------------------------------------
void dolfin::MPI::reorderComm(Mesh& mesh)
{

  message("Reordering MPI Comm");
#ifdef HAVE_MPI
  int nnodes = (int) numProcesses();
  int *index = new int[nnodes];
  short *neigh = new short[nnodes];

  memset(neigh, 0, nnodes * sizeof(short));
  for (MeshGhostIterator it(mesh.distdata(), 0); !it.end(); ++it)
    neigh[it.owner()] = 1;

  MPI_Status status;
  uint src, dest;
  short recv;

  for (uint i = 1; i < numProcesses(); i++)
  {
    src = (processNumber() - i + numProcesses()) % numProcesses();
    dest = (processNumber() + i) % numProcesses();

    MPI_Sendrecv(&neigh[dest], 1, MPI_SHORT, dest, 0, &recv, 1, MPI_SHORT, src,
                 0, DOLFIN_COMM, &status);
    neigh[src] |= recv;
  }

  int degree = 0;
  for (int i = 0; i < nnodes; i++)
    degree += neigh[i];
  dolfin_assert(degree < (nnodes - 1));
  int offset;
  MPI_Scan(&degree, &offset, 1, MPI_INT, MPI_SUM, DOLFIN_COMM);
  MPI_Allgather(&offset, 1, MPI_INT, index, 1, MPI_INT, DOLFIN_COMM);
  int nedges = 0;
  MPI_Allreduce(&degree, &nedges, 1, MPI_INT, MPI_SUM, DOLFIN_COMM);
  int *edges = new int[nedges];

  int j = (int) (processNumber() > 0 ? index[processNumber() - 1] : 0);
  for (int i = 0; i < nnodes; i++)
  {
    if (neigh[i]) edges[j++] = i;
    if (j >= index[processNumber()]) break;
  }

  int *displs = new int[nnodes];
  offset -= degree;
  MPI_Allgather(&offset, 1, MPI_INT, displs, 1, MPI_INT, DOLFIN_COMM);

  int *recvcount = new int[nnodes];
  recvcount[0] = index[0];
  for (int i = 1; i < nnodes; i++)
    recvcount[i] = index[i] - index[i - 1];

  MPI_Allgatherv(&edges[(processNumber() > 0 ? index[processNumber() - 1] : 0)],
                 recvcount[processNumber()], MPI_INT, edges, recvcount, displs,
                 MPI_INT, DOLFIN_COMM);

  MPI_Comm TMP_COMM;
  if (_dolfin_comm) MPI_Comm_dup(DOLFIN_COMM, &TMP_COMM);
  else MPI_Comm_dup(MPI_COMM_WORLD, &TMP_COMM);

  MPI_Graph_create(TMP_COMM, nnodes, index, edges, 1, &DOLFIN_COMM);

  delete[] recvcount;
  delete[] displs;
  delete[] edges;
  delete[] neigh;
  delete[] index;

  int old_rank;
  MPI_Comm_rank(TMP_COMM, &old_rank);
  MPI_Comm_rank(DOLFIN_COMM, &this_process);

  int *process_map = new int[nnodes];
  process_map[old_rank] = this_process;
  MPI_Allgather(&process_map[old_rank], 1, MPI_INT, process_map, 1, MPI_INT,
                TMP_COMM);

  for (int i = 0; i < nnodes; i++)
  {
    if (process_map[i] != i)
    {
      message("Communicator changed, rebuilding mesh structure");
      mesh.distdata().remap_ownership(process_map);
      break;
    }
  }
  delete[] process_map;

  MPI_Comm_free(&TMP_COMM);
#endif
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::seed()
{
  return this_seed;
}
//-----------------------------------------------------------------------------
bool dolfin::MPI::is_valid_rank(uint rank)
{
  return rank < static_cast<uint>(num_processes);
}
//-----------------------------------------------------------------------------
