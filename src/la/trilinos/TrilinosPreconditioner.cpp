// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_TRILINOS

#include <dolfin/la/trilinos/TrilinosPreconditioner.h>

#include <dolfin/la/PreconditionerType.h>
#include <dolfin/la/trilinos/TrilinosKrylovSolver.h>
#include <dolfin/la/trilinos/TrilinosVector.h>
#include <dolfin/la/trilinos/TrilinosMatrix.h>

namespace dolfin
{

namespace trilinos
{

//-----------------------------------------------------------------------------

Preconditioner::Preconditioner( PreconditionerType type )
  : pc_type_( type )
{
}

//-----------------------------------------------------------------------------

auto Preconditioner::init( Matrix const & P, std::string pc_str ) -> void
{
  init( P, pc_type( pc_str ) );
}

//-----------------------------------------------------------------------------

auto Preconditioner::init( Matrix const & P, PreconditionerType pc_type ) -> void
{
  pc_type_ = pc_type;

  init( P );
}

//-----------------------------------------------------------------------------

auto Preconditioner::init( Matrix const & P ) -> void
{
  Teuchos::RCP< Teuchos::ParameterList > params = Teuchos::parameterList();
  switch ( pc_type_ )
  {
    case cheb:
      name_ = "CHEB";
      // params->set( "chebyshev: max eigenvalue", (real) );
      // params->set( "chebyshev: ratio eigenvalue", (real, default 30) );
      // params->set( "chebyshev: min eigenvalue", (real) );
      // params->set( "chebyshev: degree", (int) );
      // params->set( "relaxation: sweeps", (int) );
      // params->set( "smoother: sweeps", (int) );
      // params->set( "chebyshev: eigenvalue max iterations", (int) );
      // params->set( "eigen-analysis: iterations", (int) );
      // params->set( "eigen-analysis: type", "power method" );
      // params->set( "chebyshev: assume matrix does not change", false );
      // params->set( "chebyshev: min diagonal value", (real) );
      // params->set( "chebyshev: zero starting solution", (bool) );
      break;
    case ilu:
      name_ = "ILUT";
      // params->set( "fact: ilut level-of-fill", (real) );
      // params->set( "fact: drop tolerance", (real) );
      break;
    case riluk:
      name_ = "RILUK";
      // params->set( "fact: iluk level-of-fill" (int) );
      // params->set( "fact: iluk overalloc" (real) );
      break;
    case sor:
      name_ = "RELAXATION";
      params->set( "relaxation: type", "Gauss-Seidel" );
      params->set( "relaxation: damping factor", 1.2 );
      // params->set( "relaxation: backward mode", false );
      break;
    case jacobi:
      name_ = "RELAXATION";
      params->set( "relaxation: type", "Jacobi" );
      params->set( "relaxation: damping factor", 1.0 );
      // params->set( "relaxation: sweeps", 1 );
      // params->set( "relaxation: zero starting solution", true );
      // params->set( "relaxation: fix tiny diagonal entries", false );
      // params->set( "relaxation: min diagonal value", (real) ); // only if above is true
      // params->set( "relaxation: check diagonal entries", false );
      break;
    case bjacobi:
    case default_pc:
      name_ = "BLOCK_RELAXATION";
      params->set( "relaxation: type", "Jacobi" );
      params->set( "relaxation: damping factor", 1.0 );
      // params->set( "relaxation: sweeps", 1 );
      // params->set( "relaxation: zero starting solution", true );
      // params->set( "relaxation: fix tiny diagonal entries", false );
      // params->set( "relaxation: min diagonal value", (real) ); // only if above is true
      // params->set( "relaxation: check diagonal entries", false );
      break;
    case none:
      name_ = "";
      break;
    default:
      error( "Preconditioner: unsupported trilinos preconditioner: %s",
             to_string( pc_type_ ).c_str() );
      break;
  }

  if ( pc_type_ != none )
  {
    Ifpack2::Factory pc;
    pc_ = pc.create( name_, P.mat().getConst() );
    pc_->setParameters( *params );
    pc_->initialize();
    pc_->compute();
  }
  else
  {
    pc_ = Teuchos::null;
  }

  message( 1, "Preconditioner: initialized \"%s\" preconditioner",
              to_string( pc_type_ ).c_str() );
}

//-----------------------------------------------------------------------------

auto Preconditioner::set( KrylovSolver & solver ) -> void
{
  if ( pc_ != Teuchos::null )
  {
    // FIXME we would like 'solver.problem_->setLeftPrec( pc_ )' here, but:
    // > Belos::BiCGStabSolMgr::solve: The left-preconditioned case has not yet
    // > been implemented.  Please use right preconditioning for now.  If you
    // > need to use left preconditioning, please contact the Belos developers.
    // > Left preconditioning is more interesting in BiCGStab because whether it
    // > works depends on the initial guess (e.g., an initial guess of all zeros
    // > might NOT work).
    // :-(
    solver.problem_->setRightPrec( pc_ );
  }
}

//-----------------------------------------------------------------------------

} // namespace trilinos

} // namespace dolfin

#endif // HAVE_TRILINOS
