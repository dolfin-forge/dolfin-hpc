// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SOLVER_TYPE_H
#define __DOLFIN_SOLVER_TYPE_H

#include <dolfin/log/log.h>

namespace dolfin
{

  /// List of predefined solvers

  enum SolverType
  {
    lu,            // LU factorization
    cg,            // Krylov conjugate gradient method
    pipecg,        // Pipelined conjugate gradient method
    asyncg,        // Gropp's asynchronous conjugate gradient method
    gmres,         // Krylov GMRES method
    bicgstab,      // Krylov stabilised biconjugate gradient squared method
    pipebicgstab,  // Pipelined stabilised biconjugate gradient squared method
    default_solver // Default Krylov solver
  };

  //---------------------------------------------------------------------------
  inline static SolverType krylov_method(std::string type)
  {
    if (type == "cg")
      return cg;
    else if (type == "gmres")
      return gmres;
    else if (type == "bicgstab")
      return bicgstab;
    else if (type == "pipebicgstab")
      return pipebicgstab;
    else
    {
      warning("Undefined solver type            "
              "Fallback to default Krylov method");
      return default_solver;
    }
  }

  //---------------------------------------------------------------------------
  inline static SolverType solver_type(std::string type, bool fallback = false)
  {
    if (type == "lu")
    {
      return lu;
    }
    else if (type == "cg")
    {
      return cg;
    }
    else if (type == "gmres")
    {
      return gmres;
    }
    else if (type == "bicgstab")
    {
      return bicgstab;
    }
    else if (type == "pipebicgstab")
    {
      return pipebicgstab;
    }
    else if (type == "default")
    {
      return default_solver;
    }
    else
    {
      if(fallback)
      {
        warning("Undefined solver type: "
                "Fallback to default Krylov method");
        return default_solver;
      }
    }
    error("Undefined solver type.");
    return default_solver;
  }

}

#endif
