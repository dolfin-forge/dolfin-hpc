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
  // : method( method )
  // , pc_petsc( pc )
{
}

//-----------------------------------------------------------------------------

KrylovSolver::KrylovSolver( SolverType                 type,
                            trilinos::Preconditioner & Preconditioner )
	: solver_type_{}
	, solver_{}
	, pc_type_{}
	, pc_ { nullptr }
	, problem_{ Teuchos::null }
{
}

//-----------------------------------------------------------------------------

KrylovSolver::~KrylovSolver()
{
}

//-----------------------------------------------------------------------------

auto KrylovSolver::solve( trilinos::Matrix const & A,
                          trilinos::Vector &       x,
                          trilinos::Vector const & b ) -> size_t
{
  size_t M = A.size( 0 );
  size_t N = A.size( 1 );

  // Check dimensions of A
  if ( M == 0 || N == 0 )
  {
    error( "KrylovSolver: Matrix does not have a nonzero number of rows and columns" );
  }

  // Check dimensions of A vs b
  if ( M != b.size() )
  {
    error( "KrylovSolver: Non-matching dimensions for linear system (matrix has"
           " %ld rows and right-hand side vector has %ld rows)", M, b.size() );
  }

  // Check dimensions of A vs x
  if ( x.size() != N )
  {
    error( "KrylovSolver: Non-matching dimensions for linear system (matrix has"
           " %ld columns and solution vector has %ld rows)", N, x.size() );
  }

  // Write a message
  if ( dolfin_get< bool >( "Krylov report" ) )
  {
    message( "Solving linear system of size %d x %d (Belos Krylov solver).", M, N );
  }

  // Reinitialize solver if necessary
  init();

  // Reinitialize solution vector if necessary
  x.init( b.local_size() );

  // preconditioner
  problem_->setOperator( A.mat() );

  if ( pc_ != nullptr	 )
  {
    pc_->init( A );
    pc_->set( *this );
  }

  // solve
  size_t num_iterations = DOLFIN_SIZE_T_MAX;

  // problem_->setProblem( x.vec(), b.vec() );
  solver_->setProblem( problem_ );
  // status

  // Update ghosts
  // FIXME does nothing so far
  // x.update_ghost_values();

  return num_iterations;
}

//-----------------------------------------------------------------------------

void KrylovSolver::disp() const
{
}

//-----------------------------------------------------------------------------

void KrylovSolver::init()
{


  Teuchos::RCP< Teuchos::ParameterList > dummy_params = Teuchos::parameterList();

  // FIXME
  std::string method_name( "GMRES" );

  Belos::SolverFactory<real, Vector::TPVector, TPOperator> factory;
  solver_  = factory.create( method_name, dummy_params );
  problem_ = Teuchos::rcp( new BLSProblem );
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
