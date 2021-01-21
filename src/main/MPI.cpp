// Copyright (C) 2007 Magnus Vikstrøm.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/main/MPI.h>

#include <dolfin/common/maybe_unused.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/SubSystemsManager.h>

#include <ctime>
#include <string>

namespace dolfin
{

//--- STATIC ------------------------------------------------------------------
MPI::Communicator MPI::DOLFIN_COMM_WORLD  = DOLFIN_COMM_NULL;
MPI::Communicator MPI::DOLFIN_COMM_SELF   = DOLFIN_COMM_NULL;
MPI::Communicator MPI::DOLFIN_COMM        = DOLFIN_COMM_NULL;

real MPI::time_ = 0.0;
bool MPI::init_ = false;
MPI::Context MPI::ctx_;

#if HAVE_MPI
#define DOLFIN_MPI_SUBSYSTEM_INIT \
  if (!init_) { SubSystemsManager::MPI::init(); }
#endif

#define DOLFIN_MPI_WRN_UNIMPLEMENTED \
  warning("Unimplemented without MPI support");

#define DOLFIN_MPI_ERR_UNIMPLEMENTED \
  error("Unimplemented without MPI support");

//-----------------------------------------------------------------------------
auto MPI::rank() -> uint
{
#ifdef HAVE_MPI
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.rank);
#else
  DOLFIN_MPI_WRN_UNIMPLEMENTED
  return 0;
#endif
}
//-----------------------------------------------------------------------------
auto MPI::size() -> uint
{
#ifdef HAVE_MPI
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.size);
#else
  DOLFIN_MPI_WRN_UNIMPLEMENTED
  return 1;
#endif
}
//-----------------------------------------------------------------------------
auto MPI::group_id() -> uint
{
#ifdef HAVE_MPI
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.group_idx);
#else
  DOLFIN_MPI_ERR_UNIMPLEMENTED
  return 0;
#endif
}
//-----------------------------------------------------------------------------
auto MPI::num_groups() -> uint
{
#ifdef HAVE_MPI
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.group_cnt);
#else
  DOLFIN_MPI_ERR_UNIMPLEMENTED
  return 0;
#endif
}
//-----------------------------------------------------------------------------
auto MPI::global_rank() -> uint
{
#ifdef HAVE_MPI
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.global_rank);
#else
  DOLFIN_MPI_ERR_UNIMPLEMENTED
  return 0;
#endif
}
//-----------------------------------------------------------------------------
auto MPI::global_size() -> uint
{
#ifdef HAVE_MPI
  DOLFIN_MPI_SUBSYSTEM_INIT
  return static_cast<uint>(ctx_.global_size);
#else
  DOLFIN_MPI_ERR_UNIMPLEMENTED
  return 0;
#endif
}
//-----------------------------------------------------------------------------
void MPI::startTimer()
{
#ifdef HAVE_MPI
  MPI::check_error( MPI_Barrier(MPI::DOLFIN_COMM) );
  time_ = MPI_Wtime();
#else
  DOLFIN_MPI_ERR_UNIMPLEMENTED
#endif
}
//-----------------------------------------------------------------------------
auto MPI::stopTimer() -> real
{
#ifdef HAVE_MPI
  MPI::check_error( MPI_Barrier(MPI::DOLFIN_COMM) );
  return (MPI_Wtime() - time_);
#else
  DOLFIN_MPI_ERR_UNIMPLEMENTED
  return 0.0;
#endif
}
//-----------------------------------------------------------------------------
void MPI::startTimer(real& stime)
{
#ifdef HAVE_MPI
  MPI::check_error( MPI_Barrier(MPI::DOLFIN_COMM) );
  stime = MPI_Wtime();
#else
  MAYBE_UNUSED(stime)
  DOLFIN_MPI_ERR_UNIMPLEMENTED
#endif
}
//-----------------------------------------------------------------------------
auto MPI::stopTimer(real& stime) -> real
{
#ifdef HAVE_MPI
  MPI::check_error( MPI_Barrier(MPI::DOLFIN_COMM) );
  return (MPI_Wtime() - stime);
#else
  MAYBE_UNUSED(stime)
  DOLFIN_MPI_ERR_UNIMPLEMENTED
  return 0.0;
#endif
}
//-----------------------------------------------------------------------------
void MPI::initComm(int ngroups)
{
  if (init_) { return; }

#ifdef HAVE_MPI

  // Initialize world
  MPI::check_error( MPI_Comm_dup(MPI_COMM_WORLD, &MPI::DOLFIN_COMM_WORLD) );
  MPI::check_error( MPI_Comm_rank(MPI::DOLFIN_COMM_WORLD, &ctx_.global_rank) );
  MPI::check_error( MPI_Comm_size(MPI::DOLFIN_COMM_WORLD, &ctx_.global_size) );
  MPI::check_error( MPI_Comm_dup(MPI_COMM_SELF, &MPI::DOLFIN_COMM_SELF) );
  ctx_.seed = std::time(nullptr) + ctx_.global_rank;

  // Initialize group(s)
  int const wsize = ctx_.global_size;
  int const wrank = ctx_.global_rank;
  if ((ngroups > 1) && (wsize >= ngroups))
  {
    MPI_Group wgroup;
    MPI::check_error( MPI_Comm_group(MPI::DOLFIN_COMM_WORLD, &wgroup) );

    int const p = wsize / ngroups;
    int const r = wsize % ngroups;
    int const k = wrank / (r ? p + 1 : p);
    int const ssize = (wrank < (p + 1) * r ? p + 1 : p);
    int const sroot = k * p + std::min(k, r);
    int srange[3] = { sroot, sroot + ssize - 1, 1 };

    MPI_Group sgroup;
    MPI::check_error( MPI_Group_range_incl(wgroup, 1, &srange, &sgroup) );
    MPI::check_error( MPI_Comm_create(MPI::DOLFIN_COMM_WORLD, sgroup,
                                      &MPI::DOLFIN_COMM) );
    MPI::check_error( MPI_Group_rank(sgroup, &ctx_.rank) );
    MPI::check_error( MPI_Group_size(sgroup, &ctx_.size) );
    ctx_.group_idx = k;
    ctx_.group_cnt = ngroups;
  }
  else
  {
    MPI::check_error( MPI_Comm_dup(MPI::DOLFIN_COMM_WORLD, &MPI::DOLFIN_COMM) );
    MPI::check_error( MPI_Comm_rank(MPI::DOLFIN_COMM, &ctx_.rank) );
    MPI::check_error( MPI_Comm_size(MPI::DOLFIN_COMM, &ctx_.size) );
    ctx_.group_idx = 0;
    ctx_.group_cnt = 1;
  }
#else
  MAYBE_UNUSED(ngroups)
#endif

  init_ = true;
}
//-----------------------------------------------------------------------------
void MPI::finiComm()
{
#ifdef HAVE_MPI
  MPI::check_error( MPI_Comm_free(&MPI::DOLFIN_COMM) );
  MPI::check_error( MPI_Comm_free(&MPI::DOLFIN_COMM_SELF) );
  MPI::check_error( MPI_Comm_free(&MPI::DOLFIN_COMM_WORLD) );
#else
  DOLFIN_MPI_ERR_UNIMPLEMENTED
#endif
}
//-----------------------------------------------------------------------------
auto MPI::seed() -> uint
{
  return ctx_.seed;
}
//-----------------------------------------------------------------------------
auto MPI::is_valid_rank(uint rank) -> bool
{
  return (rank < static_cast<uint>(ctx_.size));
}
//-----------------------------------------------------------------------------
auto MPI::is_root() -> bool
{
  return (ctx_.rank == 0);
}
//-----------------------------------------------------------------------------
void MPI::offset(uint local, uint& offset, Communicator& comm)
{
  // Fool-proof as the value for rank 0 is undefined according to MPI specs
  offset = 0;
#ifdef HAVE_MPI
  MPI::exscan_sum( &local, &offset, 1, comm );
#else
  MAYBE_UNUSED(local)
  MAYBE_UNUSED(comm)
#endif
}
//-----------------------------------------------------------------------------
#if defined( DEBUG )
auto MPI::check_error( int const mpi_error ) -> int
{
#if defined(HAVE_MPI)
  switch ( mpi_error )
  {
  case MPI_SUCCESS:
    break;
  case MPI_ERR_BUFFER:
   error( "Invalid buffer pointer" );
    break;
  case MPI_ERR_COUNT:
   error( "Invalid count argument" );
    break;
  case MPI_ERR_TYPE:
   error( "Invalid datatype argument" );
    break;
  case MPI_ERR_TAG:
   error( "Invalid tag argument" );
    break;
  case MPI_ERR_COMM:
   error( "Invalid communicator" );
    break;
  case MPI_ERR_RANK:
   error( "Invalid rank" );
    break;
  case MPI_ERR_REQUEST:
   error( "Invalid request (handle)" );
    break;
  case MPI_ERR_ROOT:
   error( "Invalid root" );
    break;
  case MPI_ERR_GROUP:
   error( "Invalid group" );
    break;
  case MPI_ERR_OP:
   error( "Invalid operation" );
    break;
  case MPI_ERR_TOPOLOGY:
   error( "Invalid topology" );
    break;
  case MPI_ERR_DIMS:
   error( "Invalid dimension argument" );
    break;
  case MPI_ERR_ARG:
   error( "Invalid argument of some other kind" );
    break;
  case MPI_ERR_UNKNOWN:
   error( "Unknown error" );
    break;
  case MPI_ERR_TRUNCATE:
   error( "Message truncated on receive" );
    break;
  case MPI_ERR_OTHER:
   error( "Known error not in this list" );
    break;
  case MPI_ERR_INTERN:
   error( "Internal MPI (implementation) error" );
    break;
  case MPI_ERR_IN_STATUS:
   error( "Error code is in status" );
    break;
  case MPI_ERR_PENDING:
   error( "Pending request" );
    break;
  case MPI_ERR_KEYVAL:
   error( "Invalid keyval has been passed" );
    break;
  case MPI_ERR_NO_MEM:
   error( "MPI_ALLOC_MEM failed because memory is exhausted" );
    break;
  case MPI_ERR_BASE:
   error( "Invalid base passed to MPI_FREE_MEM" );
    break;
  case MPI_ERR_INFO_KEY:
   error( "Key longer than MPI_MAX_INFO_KEY" );
    break;
  case MPI_ERR_INFO_VALUE:
   error( "Value longer than MPI_MAX_INFO_VAL" );
    break;
  case MPI_ERR_INFO_NOKEY:
   error( "Invalid key passed to MPI_INFO_DELETE" );
    break;
  case MPI_ERR_SPAWN:
   error( "Error in spawning processes" );
    break;
  case MPI_ERR_PORT:
   error( "Invalid port name passed to MPI_COMM_CONNECT" );
    break;
  case MPI_ERR_SERVICE:
   error( "Invalid service name passed to MPI_UNPUBLISH_NAME" );
    break;
  case MPI_ERR_NAME:
   error( "Invalid service name passed to MPI_LOOKUP_NAME" );
    break;
  case MPI_ERR_WIN:
   error( "Invalid win argument" );
    break;
  case MPI_ERR_SIZE:
   error( "Invalid size argument" );
    break;
  case MPI_ERR_DISP:
   error( "Invalid disp argument" );
    break;
  case MPI_ERR_INFO:
   error( "Invalid info argument" );
    break;
  case MPI_ERR_LOCKTYPE:
   error( "Invalid locktype argument" );
    break;
  case MPI_ERR_ASSERT:
   error( "Invalid assert argument" );
    break;
  case MPI_ERR_RMA_CONFLICT:
   error( "Conflicting accesses to window" );
    break;
  case MPI_ERR_RMA_SYNC:
   error( "Wrong synchronization of RMA calls" );
    break;
#if ( MPI_VERSION >  2)
  case MPI_ERR_RMA_RANGE:
   error( "Target memory is not part of the window (in the case of a window \
           created with MPI_WIN_CREATE_DYNAMIC, target memory not attached)" );
    break;
  case MPI_ERR_RMA_ATTACH:
   error( "Memory cannot be attached (e.g., because of resource exhaustion)" );
    break;
  case MPI_ERR_RMA_SHARED:
   error( "Memory cannot be shared (e.g., some process in the group of the \
           specified communicator cannot expose shared memory)" );
    break;
  case MPI_ERR_RMA_FLAVOR:
   error( "Passed window has the wrong flavor for the called function" );
    break;
#endif
  case MPI_ERR_FILE:
   error( "Invalid file handle" );
    break;
  case MPI_ERR_NOT_SAME:
   error( "Collective argument not identical on all processes, or collective \
           routines called in a different order by different processes" );
    break;
  case MPI_ERR_AMODE:
   error( "Error related to the amode passed to MPI_FILE_OPEN" );
    break;
  case MPI_ERR_UNSUPPORTED_DATAREP:
   error( "Unsupported datarep passed to MPI_FILE_SET_VIEW" );
    break;
  case MPI_ERR_UNSUPPORTED_OPERATION:
   error( "Unsupported operation, such as seeking on a file which supports \
           sequential access only" );
    break;
  case MPI_ERR_NO_SUCH_FILE:
   error( "File does not exist" );
    break;
  case MPI_ERR_FILE_EXISTS:
   error( "File exists" );
    break;
  case MPI_ERR_BAD_FILE:
   error( "Invalid file name (e.g., path name too long)" );
    break;
  case MPI_ERR_ACCESS:
   error( "Permission denied" );
    break;
  case MPI_ERR_NO_SPACE:
   error( "Not enough space" );
    break;
  case MPI_ERR_QUOTA:
   error( "Quota exceeded" );
    break;
  case MPI_ERR_READ_ONLY:
   error( "Read-only file or file system" );
    break;
  case MPI_ERR_FILE_IN_USE:
   error( "File operation could not be completed, as the file is currently \
           open by some process" );
    break;
  case MPI_ERR_DUP_DATAREP:
   error( "Conversion functions could not be registered because a data \
           representation identifier that was already defined was passed to \
           MPI_REGISTER_DATAREP" );
    break;
  case MPI_ERR_CONVERSION:
   error( "An error occurred in a user supplied data conversion function." );
    break;
  case MPI_ERR_IO:
   error( "Other I/O error" );
    break;
  case MPI_ERR_LASTCODE:
   error( "Last error code" );
   break;
  }
#endif

  return mpi_error;
}
#endif // DEBUG
#ifdef HAVE_MPI
//-----------------------------------------------------------------------------
void MPI::file_open( MPI_File & file, std::string const & filename,
                    int mode, Communicator & comm, MPI_Info info )
{
  int const error = MPI_File_open( comm, filename.c_str(), mode, info, &file );
  check_file_status( error );
  check_error( error );
}
//-----------------------------------------------------------------------------
void MPI::file_close( MPI_File & file )
{
  int const error = MPI_File_close( &file );
  check_file_status( error );
  check_error( error );
}
//-----------------------------------------------------------------------------
void MPI::check_file_status( int const mpi_error )
{
  switch ( mpi_error )
  {
  case MPI_ERR_FILE:
   error( "Invalid file handle" );
    break;
  case MPI_ERR_AMODE:
   error( "Error related to the amode passed to MPI_FILE_OPEN" );
    break;
  case MPI_ERR_UNSUPPORTED_DATAREP:
   error( "Unsupported datarep passed to MPI_FILE_SET_VIEW" );
    break;
  case MPI_ERR_UNSUPPORTED_OPERATION:
   error( "Unsupported operation, such as seeking on a file which supports \
           sequential access only" );
    break;
  case MPI_ERR_NO_SUCH_FILE:
   error( "File does not exist" );
    break;
  case MPI_ERR_FILE_EXISTS:
   error( "File exists" );
    break;
  case MPI_ERR_BAD_FILE:
   error( "Invalid file name (e.g., path name too long)" );
    break;
  case MPI_ERR_ACCESS:
   error( "Permission denied" );
    break;
  case MPI_ERR_NO_SPACE:
   error( "Not enough space" );
    break;
  case MPI_ERR_QUOTA:
   error( "Quota exceeded" );
    break;
  case MPI_ERR_READ_ONLY:
   error( "Read-only file or file system" );
    break;
  case MPI_ERR_FILE_IN_USE:
   error( "File operation could not be completed, as the file is currently \
           open by some process" );
    break;
  case MPI_ERR_IO:
   error( "Other I/O error" );
   default:
    break;
  }
}
//-----------------------------------------------------------------------------
#endif
} /* namespace dolfin */
