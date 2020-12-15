// Copyright (C) 2005 Johan Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/la/PETScKrylovSolver.h>

#include <dolfin/config/dolfin_config.h>

#include <dolfin/la/PETScKrylovMatrix.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/PETScVector.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/parameter/parameters.h>

#ifdef HAVE_PETSC

// Is this needed in PETSc 3.x?
#if PETSC_VERSION_MAJOR == 2
#include <private/pcimpl.h>
#endif

#include <petscksp.h>

namespace dolfin
{

// Monitor function
int monitor( KSP, int iteration, real rnorm, void * )
{
  message( "Iteration %d: residual = %g", iteration, rnorm );
  return 0;
}

//-----------------------------------------------------------------------------
PETScKrylovSolver::PETScKrylovSolver( SolverType method, PreconditionerType pc )
  : method( method )
  , pc_petsc( pc )
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScKrylovSolver::PETScKrylovSolver( SolverType            method,
                                      PETScPreconditioner & preconditioner )
  : method( method )
  , pc_petsc( default_pc )
  , pc_dolfin( &preconditioner )
  , ksp( nullptr )
  , M( 0 )
  , N( 0 )
  , parameters_read( false )
  , pc_set( false )
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScKrylovSolver::~PETScKrylovSolver()
{
  // Destroy solver environment.
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  if ( ksp )
    KSPDestroy( &ksp );
#else
  if ( ksp )
    KSPDestroy( ksp );
#endif
}
//-----------------------------------------------------------------------------
dolfin::uint PETScKrylovSolver::solve( const PETScMatrix & A,
                                       PETScVector &       x,
                                       const PETScVector & b )
{
  // Check dimensions
  uint M = A.size( 0 );
  uint N = A.size( 1 );
  if ( N != b.size() )
  {
    error( "Non-matching dimensions for linear system." );
  }

  // Write a message
  if ( dolfin_get< bool >( "Krylov report" ) )
  {
    message( "Solving linear system of size %d x %d (Krylov solver).", M, N );
  }

  // Reinitialize KSP solver if necessary
  init( M, N );

  // Reinitialize solution vector if necessary
  x.init( b.local_size() );

  // Read parameters if not done
  if ( !parameters_read )
    readParameters();
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 4
  KSPSetOperators( ksp, A.mat(), A.mat() );
  if ( dolfin_get< bool >( "Krylov keep PC" ) )
    KSPSetReusePreconditioner( ksp, PETSC_TRUE );

#else
  // Solve linear system
  if ( dolfin_get< bool >( "Krylov keep PC" ) )
  {
    KSPSetOperators( ksp, A.mat(), A.mat(), SAME_PRECONDITIONER );
  }
  else
  {
    KSPSetOperators( ksp, A.mat(), A.mat(), SAME_NONZERO_PATTERN );
  }
#endif

  // FIXME: Preconditioner being set here to avoid PETSc bug with Hypre.
  //        See explanation inside PETScKrylovSolver:init().
  if ( !pc_set )
  {
    setPETScPreconditioner();
    pc_set = true;
  }

  KSPSolve( ksp, b.vec(), x.vec() );

  int num_iterations = DOLFIN_INT_MAX;

  // Check if the solution converged
  KSPConvergedReason reason;
  KSPGetConvergedReason( ksp, &reason );

  std::string message = "";

  if ( reason == KSP_DIVERGED_ITS )
  {
    message = "Required more than \"its\" to reach convergence.";
  }
  else if ( reason == KSP_DIVERGED_DTOL )
  {
    message = "Residual norm increased by a factor of divtol.";
  }
  else if ( reason == KSP_DIVERGED_NANORINF )
  {
    message = "Residual norm became Not-a-number or Inf likely due to 0/0.";
  }
  else if ( reason == KSP_DIVERGED_BREAKDOWN )
  {
    message = "Generic breakdown in method.";
  }
  else if ( reason == KSP_DIVERGED_BREAKDOWN_BICG )
  {
    message = "Initial residual is orthogonal to preconditioned initial residual. "
              "Try a different preconditioner, or a different initial Level. ";
  }
  else if ( reason == KSP_DIVERGED_INDEFINITE_PC )
  {
    message = "Divergence because of indefinite preconditioner.";
  }
  else if ( reason == KSP_DIVERGED_NONSYMMETRIC )
  {
    message = "KSP_DIVERGED_NONSYMMETRIC";
  }
  else if ( reason == KSP_DIVERGED_INDEFINITE_PC )
  {
    message = "KSP_DIVERGED_INDEFINITE_PC";
  }
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 10
  else if ( reason == KSP_DIVERGED_PC_FAILED )
#else
  else if ( reason == KSP_DIVERGED_PCSETUP_FAILED )
#endif
  {
    message = "KSP_DIVERGED_PC_FAILED";
  }
  else if ( reason < 0 )
  {
    message = "Krylov solver did not converge.";
  }

  if( message != "" )
  {
    message = "KrylovSolver did not converge.\nReason: " + message;
    if ( dolfin_get< bool >( "Krylov error on nonconvergence" ) )
    {
      error( message.c_str() );
    }
    else
    {
      warning( message.c_str() );
    }
  }
  else
  {
    // Get the number of iterations
    KSPGetIterationNumber( ksp, &num_iterations );

    // Report results
    writeReport( num_iterations );
  }

  return num_iterations;
}
//-----------------------------------------------------------------------------
dolfin::uint PETScKrylovSolver::solve( const PETScKrylovMatrix & A,
                                       PETScVector &             x,
                                       const PETScVector &       b )
{
  // Check dimensions
  uint M = A.size( 0 );
  uint N = A.size( 1 );
  if ( N != b.size() )
  {
    error( "Non-matching dimensions for linear system." );
  }

  // Write a message
  if ( dolfin_get< bool >( "Krylov report" ) )
  {
    message(
      "Solving virtual linear system of size %d x %d (Krylov solver).", M, N );
  }

  // Reinitialize KSP solver if necessary
  init( M, N );

  // Reinitialize solution vector if necessary
  x.init( b.local_size() );

  // Read parameters if not done
  if ( !parameters_read )
    readParameters();

  // Don't use preconditioner that can't handle virtual (shell) matrix
  if ( !pc_dolfin )
  {
    PC pc;
    KSPGetPC( ksp, &pc );
    PCSetType( pc, PCNONE );
  }

  // Solve linear system
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 4
  KSPSetOperators( ksp, A.mat(), A.mat() );
#else
  KSPSetOperators( ksp, A.mat(), A.mat(), DIFFERENT_NONZERO_PATTERN );
#endif
  KSPSolve( ksp, b.vec(), x.vec() );

  int num_iterations = DOLFIN_INT_MAX;

  // Check if the solution converged
  KSPConvergedReason reason;
  KSPGetConvergedReason( ksp, &reason );

  std::string message = "";

  if ( reason == KSP_DIVERGED_ITS )
  {
    message = "Required more than \"its\" to reach convergence.";
  }
  else if ( reason == KSP_DIVERGED_DTOL )
  {
    message = "Residual norm increased by a factor of divtol.";
  }
  else if ( reason == KSP_DIVERGED_NANORINF )
  {
    message = "Residual norm became Not-a-number or Inf likely due to 0/0.";
  }
  else if ( reason == KSP_DIVERGED_BREAKDOWN )
  {
    message = "Generic breakdown in method.";
  }
  else if ( reason == KSP_DIVERGED_BREAKDOWN_BICG )
  {
    message = "Initial residual is orthogonal to preconditioned initial residual. "
              "Try a different preconditioner, or a different initial Level. ";
  }
  else if ( reason == KSP_DIVERGED_INDEFINITE_PC )
  {
    message = "Divergence because of indefinite preconditioner.";
  }
  else if ( reason == KSP_DIVERGED_NONSYMMETRIC )
  {
    message = "KSP_DIVERGED_NONSYMMETRIC";
  }
  else if ( reason == KSP_DIVERGED_INDEFINITE_PC )
  {
    message = "KSP_DIVERGED_INDEFINITE_PC";
  }
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 10
  else if ( reason == KSP_DIVERGED_PC_FAILED )
#else
  else if ( reason == KSP_DIVERGED_PCSETUP_FAILED )
#endif
  {
    message = "KSP_DIVERGED_PC_FAILED";
  }
  else if ( reason < 0 )
  {
    message = "Krylov solver did not converge.";
  }

  if( message != "" )
  {
    message = "KrylovSolver did not converge.\nReason: " + message;
    if ( dolfin_get< bool >( "Krylov error on nonconvergence" ) )
    {
      error( message.c_str() );
    }
    else
    {
      warning( message.c_str() );
    }
  }
  else
  {
    // Get the number of iterations
    KSPGetIterationNumber( ksp, &num_iterations );

    // Report results
    writeReport( num_iterations );
  }

  return num_iterations;
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::disp() const
{
  KSPView( ksp, PETSC_VIEWER_STDOUT_WORLD );
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::init( uint M, uint N )
{
  // Check if we need to reinitialize
  if ( ksp != nullptr && M == this->M && N == this->N )
    return;

  // Save size of system
  this->M = M;
  this->N = N;

  // Destroy old solver environment if necessary
  if ( ksp != nullptr )
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    KSPDestroy( &ksp );
#else
    KSPDestroy( ksp );
#endif

  // Set up solver environment
#ifdef DOLFIN_HAVE_MPI
  if ( MPI::size() > 1 )
    KSPCreate( MPI::DOLFIN_COMM, &ksp );
  else
    KSPCreate( PETSC_COMM_SELF, &ksp );
#else
  KSPCreate( PETSC_COMM_SELF, &ksp );
#endif
  KSPSetFromOptions( ksp );
  KSPSetInitialGuessNonzero( ksp, PETSC_TRUE );

  // Set solver
  setSolver();

  // FIXME: The preconditioner is being set in solve() due to a PETSc bug
  //        when using Hypre preconditioner. The problem can be avoided by
  //        setting the preconditioner after KSPSetOperators(). This will be
  //        fixed in PETSc, the the preconditioner can be set here again.
  // Set preconditioner
  //  setPETScPreconditioner();
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::readParameters()
{
  // Don't do anything if not initialized
  if ( !ksp )
    return;

  // Set tolerances
  KSPSetTolerances( ksp,
                    dolfin_get< real >( "Krylov relative tolerance" ),
                    dolfin_get< real >( "Krylov absolute tolerance" ),
                    dolfin_get< real >( "Krylov divergence limit" ),
                    dolfin_get< uint >( "Krylov maximum iterations" ) );

  // Set nonzero shift for preconditioner
  if ( !pc_dolfin )
  {
    PC pc;
    KSPGetPC( ksp, &pc );

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 0
    PCFactorSetShiftType( pc, MAT_SHIFT_NONZERO );
    PCFactorSetShiftAmount( pc, dolfin_get< real >( "Krylov shift nonzero" ) );
#else
    PCFactorSetShiftNonzero( pc, dolfin_get< real >( "Krylov shift nonzero" ) );
#endif
  }

  // Remember that we have read parameters
  parameters_read = true;
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::setSolver()
{
  // Don't do anything for default method
  if ( method == default_solver )
    return;

  // Set PETSc Krylov solver
  KSPType ksp_type = getType( method );
  KSPSetType( ksp, ksp_type );
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::setPETScPreconditioner()
{
  // Treat special case DOLFIN user-defined preconditioner
  if ( pc_dolfin )
  {
    PETScPreconditioner::setup( ksp, *pc_dolfin );
    return;
  }

  // Treat special case default preconditioner (do nothing)
  if ( pc_petsc == default_pc )
    return;

  // Get PETSc PC pointer
  PC pc;
  KSPGetPC( ksp, &pc );

  // Make sure options are set
  PCSetFromOptions( pc );

  // Treat special case Hypre AMG preconditioner
  if ( pc_petsc == amg )
  {
#if PETSC_HAVE_HYPRE
    // PCSetFromOptions(pc);
    // pc->ops->setfromoptions(pc);
    PCSetType( pc, PCHYPRE );
    // pc->ops->setfromoptions(pc);
    PCHYPRESetType( pc, "boomeramg" );
    PCSetFromOptions( pc );

    // pc->ops->setfromoptions(pc);
    // PCHYPRESetFromOptions(pc);
#else
    warning( "PETSc has not been compiled with the HYPRE library for   "
             "algerbraic multigrid. Default PETSc solver will be used. "
             "For performance, installation of HYPRE is recommended.   "
             "See the DOLFIN user manual for more information." );
#endif
    return;
  }

  // Set preconditioner
  PCSetType( pc, PETScPreconditioner::getType( pc_petsc ) );
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::writeReport( int num_iterations )
{
  // Check if we should write the report
  if ( dolfin_get< bool >( "Krylov report" ) )
  {
#if PETSC_VERSION_MAJOR > 2
#if PETSC_VERSION_MINOR > 3
    KSPType ksp_type;
    PCType  pc_type;
#else
    const KSPType ksp_type;
    const PCType  pc_type;
#endif
#else
    KSPType ksp_type;
    PCType  pc_type;
#endif

    // Get name of solver
    KSPGetType( ksp, &ksp_type );

    // Get name of preconditioner
    PC pc;
    KSPGetPC( ksp, &pc );
    PCGetType( pc, &pc_type );

    // Report number of iterations and solver type
    message( "Krylov solver (%s, %s) converged in %d iterations.",
             ksp_type,
             pc_type,
             num_iterations );
  }
}
//-----------------------------------------------------------------------------
KSPType PETScKrylovSolver::getType( SolverType method ) const
{
  switch ( method )
  {
    case bicgstab:
      return KSPBCGS;
    case cg:
      return KSPCG;
    case default_solver:
      return "default";
    case gmres:
      return KSPGMRES;
    default:
      warning( "Requested Krylov method unknown. Using GMRES." );
      return KSPGMRES;
  }
}
//-----------------------------------------------------------------------------

}

#endif
