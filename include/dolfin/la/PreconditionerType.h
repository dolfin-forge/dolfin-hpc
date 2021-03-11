// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRECONDITIONER_TYPE_H
#define __DOLFIN_PRECONDITIONER_TYPE_H

#include <dolfin/log/log.h>

namespace dolfin
{

//---------------------------------------------------------------------------

/// List of predefined preconditioners

enum class PreconditionerType
{
  none,      // No preconditioning
  jacobi,    // Jacobi
  bjacobi,   // Block Jacobi
  sor,       // SOR (successive over relaxation)
  ilu,       // Incomplete LU factorization
  dilu,      // diagonal Incomplete LU factorization
  icc,       // Incomplete Cholesky factorization
  amg,       // Algebraic multigrid (through Hypre when available)
  default_pc // Default choice of preconditioner
};

//---------------------------------------------------------------------------

inline static auto pc_type( std::string type ) -> PreconditionerType
{
  if ( type == "none" )
    return PreconditionerType::none;
  else if ( type == "bjacobi" )
    return PreconditionerType::bjacobi;
  else if ( type == "sor" )
    return PreconditionerType::sor;
  else if ( type == "ilu" )
    return PreconditionerType::ilu;
  else if ( type == "dilu" )
    return PreconditionerType::dilu;
  else if ( type == "amg" )
    return PreconditionerType::amg;
  else
  {
    warning( "Undefined preconditioner          "
             "Fallback to default preconditioner" );
    return PreconditionerType::default_pc;
  }
}

//---------------------------------------------------------------------------

} // namespace dolfin

#endif
