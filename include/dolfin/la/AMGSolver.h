// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_AMG_SOLVER_H
#define __DOLFIN_AMG_SOLVER_H

#include <dolfin/common/Timer.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/MultigridScheme.h>
#include <dolfin/la/janpack/JANPACKAMGSolver.h>
#include <dolfin/la/janpack/JANPACKMat.h>
#include <dolfin/la/janpack/JANPACKVec.h>

namespace dolfin
{

/// This class defines an interface for a AMG solver.

class AMGSolver
{
public:
  /// Create Krylov solver
  AMGSolver( MultigridScheme     scheme_type     = default_scheme,
             MultigridSmoother   smoother_type   = default_smoother,
             MultigridCoarsening coarsening_type = default_coarsening )
    : scheme_type( scheme_type )
    , smoother_type( smoother_type )
    , coarsening_type( coarsening_type )
    , janpack_solver( 0 )
  {
    MAYBE_UNUSED( scheme_type );
    MAYBE_UNUSED( smoother_type );
    MAYBE_UNUSED( coarsening_type );
  }

  /// Destructor
  ~AMGSolver()
  {
    delete janpack_solver;
  }

  /// Solve linear system Ax = b
  size_t
    solve( const GenericMatrix & A, GenericVector & x, const GenericVector & b )
  {
    Timer timer( "AMG solver" );

#if defined( HAVE_JANPACK ) && !defined( HAVE_JANPACK_MPI )
    if ( A.has_type< JANPACKMat >() )
    {
      if ( !janpack_solver )
      {
        janpack_solver =
          new JANPACKAMGSolver( scheme_type, smoother_type, coarsening_type );
      }
      return janpack_solver->solve( A.down_cast< JANPACKMat >(),
                                    x.down_cast< JANPACKVec >(),
                                    b.down_cast< JANPACKVec >() );
    }
#else
    MAYBE_UNUSED( A );
    MAYBE_UNUSED( x );
    MAYBE_UNUSED( b );
    MAYBE_UNUSED( scheme_type );
    MAYBE_UNUSED( smoother_type );
    MAYBE_UNUSED( coarsening_type );
#endif
    error( "No AMG solver for given backend" );
    return 0;
  }

private:
  // Multigrid Scheme
  MultigridScheme scheme_type;

  // Multigrid Smoother
  MultigridSmoother smoother_type;

  // Multigrid Coarsening scheme
  MultigridCoarsening coarsening_type;

#if defined( HAVE_JANPACK ) && !defined( HAVE_JANPACK_MPI )
  JANPACKAMGSolver * janpack_solver;
#else
  int * janpack_solver;
#endif
};
}

#endif
