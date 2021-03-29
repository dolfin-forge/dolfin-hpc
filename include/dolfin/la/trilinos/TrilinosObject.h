// Copyright (C) 2020 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TRILINOS_OBJECT_H
#define __DOLFIN_TRILINOS_OBJECT_H

#include <dolfin/main/SubSystemsManager.h>
#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#ifdef DOLFIN_HAVE_MPI
#include <Teuchos_Comm.hpp>
#include <Tpetra_Core.hpp>
#include <dolfin/main/MPI.h>
#endif

namespace dolfin
{

namespace trilinos
{

/// This class calls SubSystemsManger to initialise Trilinos.
///
/// All Trilinos objects must be derived from this class.

using MPIComm = Teuchos::RCP< Teuchos::MpiComm< int > const >;

class Object
{
public:
  Object()
#if DOLFIN_HAVE_MPI
  	// FIXME use the correct communicator here
  	: comm_( new Teuchos::MpiComm< int >( Teuchos::MpiComm< int >( MPI::DOLFIN_COMM ) ) )
#endif
  {
    SubSystemsManager::Trilinos::initialize();
  }

protected:
  MPIComm comm_;
};

} // end namespace trilinos

} // end namespace dolfin

#endif // HAVE_TRILINOS

#endif // __DOLFIN_TRILINOS_OBJECT_H
