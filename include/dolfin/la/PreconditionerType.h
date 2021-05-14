// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRECONDITIONER_TYPE_H
#define __DOLFIN_PRECONDITIONER_TYPE_H

#include <dolfin/log/log.h>

namespace dolfin
{

//---------------------------------------------------------------------------

/// List of predefined preconditioners

enum PreconditionerType
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
  riluk,     // Relaxed ILU with level k fill (Trilinos/Ifpack2)
  default_pc // Default choice of preconditioner
};

//---------------------------------------------------------------------------

inline static auto pc_type( std::string type ) -> PreconditionerType
{
  if ( type == "none" )
    return none;
  else if ( type == "bjacobi" )
    return bjacobi;
  else if ( type == "sor" )
    return sor;
  else if ( type == "ilu" )
    return ilu;
  else if ( type == "dilu" )
    return dilu;
  else if ( type == "amg" )
    return amg;
  else if ( type == "cheb" )
    return cheb;
  else if ( type == "riluk" )
    return riluk;
  else
  {
    warning( "Undefined preconditioner          "
             "Fallback to default preconditioner" );
    return default_pc;
  }
}


//---------------------------------------------------------------------------

inline static auto to_string( PreconditionerType type ) -> std::string
{
  std::string name = "";

  switch( type )
  {
    case none:
      name = "none";
      break;
    case jacobi:
      name = "jacobi";
      break;
    case bjacobi:
      name = "bjacobi";
      break;
    case sor:
      name = "sor";
      break;
    case ilu:
      name = "ilu";
      break;
    case dilu:
      name = "dilu";
      break;
    case icc:
      name = "icc";
      break;
    case amg:
      name = "amg";
      break;
    case cheb:
      name = "cheb";
      break;
    case riluk:
      name = "riluk";
      break;
    case default_pc:
      name = "default_pc";
      break;
  }

  return name;
}

//---------------------------------------------------------------------------

} // namespace dolfin

#endif
