// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __MULTIGRID_SCHEME_H
#define __MULTIGRID_SCHEME_H

namespace dolfin
{

  /// List of predefined multigrid schemes

  enum MultigridScheme
  {
    v,                  // V-cycle scheme
    w,                  // W-cycle scheme
    fmv,                // Full Multigrid V-cycle
    default_scheme      // Default multigrid scheme
  };

  enum MultigridSmoother
  {
    mg_jacobi,          // Jacobi
    mg_gauss_seidel,    // Gauss-Seidel 
    mg_sor,             // SOR (successive over relaxation)
    default_smoother    // Default smoother
  };

}

#endif
