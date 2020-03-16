// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/la/LinearSolver.h>

#include <dolfin/la/KrylovSolver.h>
#include <dolfin/la/LUSolver.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
LinearSolver::LinearSolver( SolverType solver_type, PreconditionerType pc_type )
  : lu_solver( NULL )
  , krylov_solver( NULL )
{
  if ( solver_type == lu )
  {
    lu_solver = new LUSolver();
  }
  else
  {
    krylov_solver = new KrylovSolver( solver_type, pc_type );
  }
}
//-----------------------------------------------------------------------------
LinearSolver::~LinearSolver()
{
  delete lu_solver;
  delete krylov_solver;
}
//-----------------------------------------------------------------------------
uint LinearSolver::solve( GenericMatrix const & A,
                          GenericVector &       x,
                          GenericVector const & b )
{
  dolfin_assert( lu_solver || krylov_solver );
  if ( lu_solver )
  {
    return lu_solver->solve( A, x, b );
  }
  else
  {
    return krylov_solver->solve( A, x, b );
  }
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
