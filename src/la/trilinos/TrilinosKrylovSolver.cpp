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

#include <BelosBiCGStabSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>

namespace dolfin
{

namespace trilinos
{

//-----------------------------------------------------------------------------

KrylovSolver::KrylovSolver( SolverType method, PreconditionerType pc )
  : solver_type_( method )
  , pc_type_( pc )
{
  init();
}

//-----------------------------------------------------------------------------

KrylovSolver::KrylovSolver( SolverType                 type,
                            trilinos::Preconditioner & pc )
  : solver_type_( type )
  , solver_ {}
  , pc_type_ {}
  , pc_( &pc )
{
  init();
}

//-----------------------------------------------------------------------------

KrylovSolver::~KrylovSolver()
{
  delete pc_;
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
    message(
      "Solving linear system of size %d x %d (Belos Krylov solver).", M, N );
  }

  // Reinitialize solver if necessary
  init();

  // Reinitialize solution vector if necessary
  x.init( b.local_size() );

  // set Matrix
  problem_->setOperator( A.mat() );

  // preconditioner
  if ( pc_ != nullptr )
  {
    pc_->init( A, pc_type_ );
    pc_->set( *this );
  }

  // setup the solver and problem
  problem_->setProblem( x.vec(), b.vec() );
  solver_->setProblem( problem_ );

  // solve
  Belos::ReturnType result = solver_->solve();
  size_t const      iters  = solver_->getNumIters();

  writeReport( result, iters );

  // Update ghosts
  // FIXME does nothing so far
  // x.update_ghost_values();

  return iters;
}

//-----------------------------------------------------------------------------

void KrylovSolver::disp() const
{
  section( "Trilinos Krylov Solver" );
  prm( "Solver type", to_string( solver_type_ ) );
  prm( "PC type",     to_string( pc_type_ ) );
  solver_->description();
  end();
}

//-----------------------------------------------------------------------------

void KrylovSolver::init()
{
  // initialize the preconditioner
  if ( pc_ == nullptr and pc_type_ != PreconditionerType::none )
  {
    pc_ = new Preconditioner;
  }

  // initialize problem
  if ( problem_ == Teuchos::null )
  {
    problem_ = Teuchos::rcp( new BLSProblem );
  }

  // initialize solver
  if ( solver_ == Teuchos::null )
  {
    Teuchos::RCP< Teuchos::ParameterList > params = Teuchos::parameterList();

    // maximum number of iterations the solver is allowed to perform (int, Default: 1000)
    params->set( "Maximum Iterations", dolfin_get< int >( "Krylov maximum iterations" ) );

    // the level that residual norms must reach to decide convergence (real, Default: 1e-8)
    params->set( "Convergence Tolerance", dolfin_get< real >( "Krylov relative tolerance" ) );

    // output verbosity (MsgType, Default: Belos::Errors)
    if ( verbose() <= 0 )
    {
      params->set( "Verbosity", Belos::Errors );
    }
    else
    {
      params->set( "Verbosity",   Belos::Warnings
                                | Belos::IterationDetails
                                | Belos::StatusTestDetails
                                | Belos::TimingDetails
                                | Belos::FinalSummary );
    }

    // style of output (OutputType, Default: Belos::General)
    params->set( "Output Style", Belos::General );

    switch( solver_type_ )
    {
      case cg:
      {
        // params->set( "Block Size", 1 );
        // params->set( "Adaptive Block Size", true );
        // params->set( "Orthogonalization", "ICGS" );
        // params->set( "Orthogonalization Constant", -1 );
        // params->set( "Show Maximum Residual Norm Only", false );
        // params->set( "Implicit Residual Scaling", false );
        // params->set( "Assert Positive Definiteness", true );
        using manager = Belos::BlockCGSolMgr< real, Vector::TPVector, TPOperator >;
        solver_ = Teuchos::rcp( new manager( problem_, params ) );
        break;
      }
      case bicgstab:
      case pipebicgstab:
      {
        using manager = Belos::BiCGStabSolMgr< real, Vector::TPVector, TPOperator >;
        solver_ = Teuchos::rcp( new manager( problem_, params ) );
        break;
      }
      case gmres:
      case default_solver:
      {
        // The maximum number of restarts. This does not include the first
        // "Num Blocks" iterations (before the first restart) (int, Default: 20)
        params->set( "Maximum Restarts", dolfin_get< int >( "Krylov GMRES restart" ) );

        // The restart length. The number of vectors (or blocks, in the case of
        // multiple right-hand sides) allocated for the Krylov basis (int, Default: 300)
        // params->set( "Num Blocks"               , 300 );

        // How to scale the implicit residual norm. The default is the norm of
        // the preconditioned initial residual (std::string)
        // params->set( "Implicit Residual Scaling", );

        // How to scale the explicit residual norm. The default is the norm of
        // the (unpreconditioned) initial residual (std::string)
        // params->set( "Explicit Residual Scaling", );

        // When solving with multiple right-hand sides: the number of right-hand
        // sides that must have converged to the given tolerance, before the
        // solver will consider all the systems converged. If -1, then the solver
        // will require that all the right-hand sides have converged before
        // declaring all the systems converged. This must be no bigger than
        // the "Block Size" parameter (int)
        // params->set( "Deflation Quorum"         , );

        // The desired orthogonalization method. Currently accepted values are
        // "DGKS", "ICGS", "IMGS", and optionally "TSQR". (std::string)
        // params->set( "Orthogonalization"        , );
        using manager = Belos::BlockGmresSolMgr< real, Vector::TPVector, TPOperator >;
        solver_ = Teuchos::rcp( new manager( problem_, params ) );
        break;
      }
      default:
      {
        error( "trilinos::KrylovSolver could not initialize solver: %s",
               to_string( solver_type_ ).c_str() );
      }
    }
  }
}

//-----------------------------------------------------------------------------

void KrylovSolver::writeReport( Belos::ReturnType result, size_t num_iterations )
{
  std::string msg = "Krylov Solver";
  msg += " (" + to_string( solver_type_ );
  msg += ", " + to_string( pc_type_ ) + ")";

  if ( result == Belos::Converged )
  {
    message( "%s converged in %d iterations",
             msg.c_str(), num_iterations );
  }
  else
  {
    if ( dolfin_get< bool >( "Krylov error on nonconvergence" ) )
    {
      error( "%s did not converge (%d iterations).",
             msg.c_str(), num_iterations );
    }
    else
    {
      warning( "%s did not converge (%d iterations).",
               msg.c_str(), num_iterations );
    }
  }
}

//-----------------------------------------------------------------------------

} // namespace trilinos

} // namespace dolfin

#endif // HAVE_TRILINOS
