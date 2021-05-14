// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/main/PE.h>

#include <dolfin/main/MPI.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
auto PE::rank() -> size_t
{
#if DOLFIN_HAVE_MPI
  if ( SubSystemsManager::active( SubSystemsManager::mpi ) )
    return MPI::rank();
#endif
  return 0;
}

//-----------------------------------------------------------------------------
auto PE::size() -> size_t
{
#if DOLFIN_HAVE_MPI
  if ( SubSystemsManager::active( SubSystemsManager::mpi ) )
    return MPI::size();
#endif
  return 1;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
