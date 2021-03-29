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
  cheb,      // Chebyshev Polynomial preconditioner (Trilinos/Ifpack2)
  ilut,      // Incomplete LU factorization with threshold (Trilinos/Ifpack2)
  riluk,     // Relaxed ILU with level k fill (Trilinos/Ifpack2)
  relax,     // Jacobi type relaxation (Trilinos/Ifpack2)
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
  else if ( type == "cheb" )
    return PreconditionerType::cheb;
  else if ( type == "ilut" )
    return PreconditionerType::ilut;
  else if ( type == "riluk" )
    return PreconditionerType::riluk;
  else if ( type == "relax" )
    return PreconditionerType::relax;
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
