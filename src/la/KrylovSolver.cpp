// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/la/KrylovSolver.h>

#include <dolfin/common/Timer.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/janpack/JANPACKMat.h>
#include <dolfin/la/janpack/JANPACKVec.h>
#include <dolfin/la/petsc/PETScMatrix.h>
#include <dolfin/la/petsc/PETScVector.h>
#include <dolfin/la/trilinos/TrilinosMatrix.h>
#include <dolfin/la/trilinos/TrilinosVector.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

KrylovSolver::KrylovSolver( SolverType solver_type, PreconditionerType pc_type )
  : solver_type( solver_type )
  , pc_type( pc_type )
{
}

//-----------------------------------------------------------------------------

KrylovSolver::~KrylovSolver()
{
  if ( janpack_solver != nullptr )
    delete janpack_solver;

  if ( petsc_solver != nullptr )
    delete petsc_solver;

  if ( trilinos_solver != nullptr )
    delete trilinos_solver;
}

//-----------------------------------------------------------------------------

auto KrylovSolver::solve( GenericMatrix const & A,
                          GenericVector &       x,
                          const GenericVector & b ) -> size_t
{
  Timer timer( "Krylov solver" );

#ifdef HAVE_PETSC
  if ( A.has_type< PETScMatrix >() )
  {
    if ( petsc_solver == nullptr )
    {
      petsc_solver = new PETScKrylovSolver( solver_type, pc_type );
    }

    return petsc_solver->solve( A.down_cast< PETScMatrix >(),
                                x.down_cast< PETScVector >(),
                                b.down_cast< PETScVector >() );
  }
#endif

#ifdef HAVE_TRILINOS
  if ( A.has_type< trilinos::Matrix >() )
  {
    if ( trilinos_solver == nullptr )
    {
      trilinos_solver = new trilinos::KrylovSolver( solver_type, pc_type );
    }

    return trilinos_solver->solve( A.down_cast< trilinos::Matrix >(),
                                   x.down_cast< trilinos::Vector >(),
                                   b.down_cast< trilinos::Vector >() );
  }
#endif

#ifdef HAVE_JANPACK
  if ( A.has_type< JANPACKMat >() )
  {
    if ( janpack_solver == nullptr )
    {
      janpack_solver = new JANPACKKrylovSolver( solver_type, pc_type );
    }

    return janpack_solver->solve( A.down_cast< JANPACKMat >(),
                                  x.down_cast< JANPACKVec >(),
                                  b.down_cast< JANPACKVec >() );
  }
#endif
  error( "No default Krylov solver for given backend" );
  return 0;
}

//-----------------------------------------------------------------------------

} // namespace dolfin
