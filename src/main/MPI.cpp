// Copyright (C) 2007 Magnus Vikstrøm.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007, 2008.
// Modified by Anders Logg, 2007.
// Modified by Niclas Jansson, 2009-2010.
//
// First added:  2007-11-30
// Last changed: 2010-01-13

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/SubSystemsManager.h>
#include <dolfin/main/MPI.h>
#include <cstring>

//-----------------------------------------------------------------------------
#ifdef HAVE_MPI
dolfin::uint dolfin::MPI::processNumber()
{
  if(!_this_process) {
    SubSystemsManager::initMPI();
    initComm();
    MPI_Comm_rank(MPI::DOLFIN_COMM, &this_process);
    _this_process = true;
  }

  return static_cast<uint>(this_process);
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::numProcesses()
{
  if(!_num_processes) {
    SubSystemsManager::initMPI();
    initComm();
    MPI_Comm_size(MPI::DOLFIN_COMM, &num_processes);
    _num_processes = true;
  }

  return static_cast<uint>(num_processes);
}
//-----------------------------------------------------------------------------
void dolfin::MPI::startTimer()
{
  MPI_Barrier(MPI::DOLFIN_COMM);
  start_time = MPI_Wtime();
}
//-----------------------------------------------------------------------------
dolfin::real dolfin::MPI::stopTimer()
{
  MPI_Barrier(MPI::DOLFIN_COMM);
  return (MPI_Wtime() - start_time);
}
//-----------------------------------------------------------------------------
void dolfin::MPI::startTimer(real& stime)
{
  MPI_Barrier(MPI::DOLFIN_COMM);
  stime = MPI_Wtime();
}
//-----------------------------------------------------------------------------
dolfin::real dolfin::MPI::stopTimer(real& stime)
{
  MPI_Barrier(MPI::DOLFIN_COMM);
  return (MPI_Wtime() - stime);
}
//-----------------------------------------------------------------------------
void dolfin::MPI::initComm()
{
  if(_dolfin_comm) 
    return;
  
  MPI_Comm_dup(MPI_COMM_WORLD, &DOLFIN_COMM);
  _dolfin_comm = true;
}
//-----------------------------------------------------------------------------
void dolfin::MPI::reorderComm(Mesh& mesh)
{

  message("Reordering MPI Comm");

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
    dest = (processNumber() + i ) % numProcesses();
    
    MPI_Sendrecv(&neigh[dest], 1, MPI_SHORT, dest, 0, 
		 &recv, 1, MPI_SHORT, src, 0,  DOLFIN_COMM, &status);
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
  for(int i = 0; i < nnodes; i++)
  {
    if (neigh[i])
      edges[j++] = i ;
    if( j >= index[processNumber()] ) break;
  }
  
  int *displs =  new int[nnodes];
  offset -= degree;
  MPI_Allgather(&offset, 1, MPI_INT, displs, 1, MPI_INT, DOLFIN_COMM);
  
  int *recvcount = new int[nnodes];  
  recvcount[0] = index[0];
  for (int i = 1; i < nnodes; i++)
    recvcount[i] = index[i] - index[i - 1];

  MPI_Allgatherv(&edges[(processNumber() > 0 ? index[processNumber() - 1] : 0)],
		 recvcount[processNumber()], MPI_INT, 
		 edges, recvcount, displs, MPI_INT, DOLFIN_COMM);

  MPI_Comm TMP_COMM;
  if(_dolfin_comm)
    MPI_Comm_dup(DOLFIN_COMM, &TMP_COMM);
  else
    MPI_Comm_dup(MPI_COMM_WORLD, &TMP_COMM);
  

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
  MPI_Allgather(&process_map[old_rank], 1, MPI_INT,process_map , 1, MPI_INT, TMP_COMM);
  
  for (int i = 0; i < nnodes; i++)
  {
    if (process_map[i] != i)
    {
      message("Communicator changed, rebuilding mesh structure");
      mesh.distdata().remap_owner(process_map);
      break;
    }
  }
  delete[] process_map;

  MPI_Comm_free(&TMP_COMM);



}
//-----------------------------------------------------------------------------
dolfin::real dolfin::MPI::start_time = 0.0;
int dolfin::MPI::this_process, dolfin::MPI::num_processes;
bool dolfin::MPI::_this_process = false;
bool dolfin::MPI::_num_processes = false;
bool dolfin::MPI::_dolfin_comm = false;
MPI_Comm dolfin::MPI::DOLFIN_COMM;
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

#endif
