// Copyright (C) 2007 Magnus Vikstrøm.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2007, 2008.
// Modified by Anders Logg, 2007.
//
// First added:  2007-11-30
// Last changed: 2008-01-07

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include "MPI.h"
#include "SubSystemsManager.h"

#ifdef HAS_MPI
  #include <mpi.h>
#endif

//-----------------------------------------------------------------------------
#ifdef HAS_MPI
dolfin::uint dolfin::MPI::processNumber()
{
  if(!_this_process) {
    SubSystemsManager::initMPI();
    MPI_Comm_rank(MPI_COMM_WORLD, &this_process);
    _this_process = true;
  }

  return static_cast<uint>(this_process);
}
//-----------------------------------------------------------------------------
dolfin::uint dolfin::MPI::numProcesses()
{
  if(!_num_processes) {
    SubSystemsManager::initMPI();
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    _num_processes = true;
  }

  return static_cast<uint>(num_processes);
}
//-----------------------------------------------------------------------------
bool dolfin::MPI::broadcast()
{
  // Always broadcast from processor number 0
  return numProcesses() > 1 && processNumber() == 0;
}
//-----------------------------------------------------------------------------
bool dolfin::MPI::receive()
{
  // Always receive on processors with numbers > 0
  return numProcesses() > 1 && processNumber() > 0;
}
//-----------------------------------------------------------------------------
void dolfin::MPI::startTimer()
{
  MPI_Barrier(MPI_COMM_WORLD);
  start_time = MPI_Wtime();
}
//-----------------------------------------------------------------------------
dolfin::real dolfin::MPI::stopTimer()
{
  MPI_Barrier(MPI_COMM_WORLD);
  return (MPI_Wtime() - start_time);
}
//-----------------------------------------------------------------------------
void dolfin::MPI::startTimer(real& stime)
{
  MPI_Barrier(MPI_COMM_WORLD);
  stime = MPI_Wtime();
}
//-----------------------------------------------------------------------------
dolfin::real dolfin::MPI::stopTimer(real& stime)
{
  MPI_Barrier(MPI_COMM_WORLD);
  return (MPI_Wtime() - stime);
}
//-----------------------------------------------------------------------------
dolfin::real dolfin::MPI::start_time = 0.0;
bool dolfin::MPI::_this_process = false;
bool dolfin::MPI::_num_processes = false;
int dolfin::MPI::this_process, dolfin::MPI::num_processes;
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
bool dolfin::MPI::broadcast()
{
  return false;
}
//-----------------------------------------------------------------------------
bool dolfin::MPI::receive()
{
  return false;
}
//-----------------------------------------------------------------------------

#endif
