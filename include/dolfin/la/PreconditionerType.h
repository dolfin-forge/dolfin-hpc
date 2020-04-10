// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRECONDITIONER_TYPE_H
#define __DOLFIN_PRECONDITIONER_TYPE_H

#include <dolfin/log/log.h>

namespace dolfin
{

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
    default_pc // Default choice of preconditioner
  };

  //---------------------------------------------------------------------------
  inline static PreconditionerType pc_type(std::string type)
  {
    if (type == "none")
      return none;
    else if (type == "bjacobi")
      return bjacobi;
    else if (type == "sor")
      return sor;
    else if (type == "ilu")
      return ilu;
    else if (type == "dilu")
      return dilu;
    else if (type == "amg")
      return amg;
    else
    {
      warning("Undefined preconditioner          "
              "Fallback to default preconditioner");
      return default_pc;
    }
  }

}

#endif
