// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosKrylovSolver.h>

#include <dolfin/la/trilinos/TrilinosMatrix.h>
#include <dolfin/la/trilinos/TrilinosVector.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/parameter/parameters.h>

namespace dolfin
{

namespace trilinos
{

//-----------------------------------------------------------------------------

KrylovSolver::KrylovSolver( SolverType method, PreconditionerType pc )
  : method( method )
  , pc_petsc( pc )
{
}

//-----------------------------------------------------------------------------

KrylovSolver::~KrylovSolver()
{
}

//-----------------------------------------------------------------------------

auto KrylovSolver::solve( trilinos::Matrix const & A,
                          trilinos::Vector &       x,
                          trilinos::Vector const & b ) -> dolfin::size_t
{
}

//-----------------------------------------------------------------------------

void KrylovSolver::disp() const
{
}

//-----------------------------------------------------------------------------

void KrylovSolver::init( size_t M, size_t N )
{
}

//-----------------------------------------------------------------------------

void KrylovSolver::readParameters()
{
}

//-----------------------------------------------------------------------------

void KrylovSolver::setSolver()
{
}

//-----------------------------------------------------------------------------

void KrylovSolver::setPETScPreconditioner()
{
}

//-----------------------------------------------------------------------------

void KrylovSolver::writeReport( int num_iterations )
{
}

//-----------------------------------------------------------------------------

} // namespace trilinos

} // namespace dolfin

#endif // HAVE_TRILINOS
