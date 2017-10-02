// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/main/PE.h>

#include <dolfin/main/MPI.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
uint PE::rank()
{
#if HAVE_MPI
  if (SubSystemsManager::active(SubSystemsManager::mpi)) return MPI::rank();
#endif
  return 0;
}

//-----------------------------------------------------------------------------
uint PE::size()
{
#if HAVE_MPI
  if (SubSystemsManager::active(SubSystemsManager::mpi)) return MPI::size();
#endif
  return 1;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
